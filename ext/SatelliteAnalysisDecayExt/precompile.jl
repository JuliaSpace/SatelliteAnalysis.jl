## Description #############################################################################
#
# Precompilation workload for the decay analysis extension.
#
############################################################################################

@setup_workload begin
    # Minimal ICGEM model with the EGM2008 coefficients up to degree 8, avoiding any
    # network access during the precompilation. The null degree-1 coefficients are
    # explicitly included since the parser leaves missing entries uninitialized.
    icgem = """
        product_type                gravity_field
        modelname                   EGM2008
        earth_gravity_constant      0.3986004415E+15
        radius                      0.63781363E+07
        max_degree                  8
        errors                      no
        norm                        fully_normalized
        tide_system                 tide_free
        end_of_head ======================================================================
        gfc 0 0 1.0d0 0.0d0
        gfc 1 0 0.0 0.0
        gfc 1 1 0.0 0.0
        gfc 2 0 -0.484165143790815e-03 0.000000000000000e+00
        gfc 2 1 -0.206615509074176e-09 0.138441389137979e-08
        gfc 2 2 0.243938357328313e-05 -0.140027370385934e-05
        gfc 3 0 0.957161207093473e-06 0.000000000000000e+00
        gfc 3 1 0.203046201047864e-05 0.248200415856872e-06
        gfc 3 2 0.904787894809528e-06 -0.619005475177618e-06
        gfc 3 3 0.721321757121568e-06 0.141434926192941e-05
        gfc 4 0 0.539965866638991e-06 0.000000000000000e+00
        gfc 4 1 -0.536157389388867e-06 -0.473567346518086e-06
        gfc 4 2 0.350501623962649e-06 0.662480026275829e-06
        gfc 4 3 0.990856766672321e-06 -0.200956723567452e-06
        gfc 4 4 -0.188519633023033e-06 0.308803882149194e-06
        gfc 5 0 0.686702913736681e-07 0.000000000000000e+00
        gfc 5 1 -0.629211923042529e-07 -0.943698073395769e-07
        gfc 5 2 0.652078043176164e-06 -0.323353192540522e-06
        gfc 5 3 -0.451847152328843e-06 -0.214955408306046e-06
        gfc 5 4 -0.295328761175629e-06 0.498070550102351e-07
        gfc 5 5 0.174811795496002e-06 -0.669379935180165e-06
        gfc 6 0 -0.149953927978527e-06 0.000000000000000e+00
        gfc 6 1 -0.759210081892527e-07 0.265122593213647e-07
        gfc 6 2 0.486488924604690e-07 -0.373789324523752e-06
        gfc 6 3 0.572451611175653e-07 0.895201130010730e-08
        gfc 6 4 -0.860237937191611e-07 -0.471425573429095e-06
        gfc 6 5 -0.267166423703038e-06 -0.536493151500206e-06
        gfc 6 6 0.947068749756882e-08 -0.237382353351005e-06
        gfc 7 0 0.905120844521618e-07 0.000000000000000e+00
        gfc 7 1 0.280887555776673e-06 0.951259362869275e-07
        gfc 7 2 0.330407993702235e-06 0.929969290624092e-07
        gfc 7 3 0.250458409225729e-06 -0.217118287729610e-06
        gfc 7 4 -0.274993935591631e-06 -0.124058403514343e-06
        gfc 7 5 0.164773255934658e-08 0.179281782751438e-07
        gfc 7 6 -0.358798423464889e-06 0.151798257443669e-06
        gfc 7 7 0.150746472872675e-08 0.241068767286303e-07
        gfc 8 0 0.494756003005199e-07 0.000000000000000e+00
        gfc 8 1 0.231607991248329e-07 0.588974540927606e-07
        gfc 8 2 0.800143604736599e-07 0.652805043667369e-07
        gfc 8 3 -0.193745381715290e-07 -0.859639339125694e-07
        gfc 8 4 -0.244360480007096e-06 0.698072508472777e-07
        gfc 8 5 -0.257011477267991e-07 0.892034891745881e-07
        gfc 8 6 -0.659648680031408e-07 0.308946730783065e-06
        gfc 8 7 0.672569751771483e-07 0.748686063738231e-07
        gfc 8 8 -0.124022771917136e-06 0.120551889384997e-06
        """

    @compile_workload begin
        path = tempname() * ".gfc"
        write(path, icgem)

        gm = GravityModels.load(IcgemFile, path)

        jd₀ = date_to_jd(2024, 1, 1)

        orb = KeplerianElements(
            jd₀,
            EARTH_EQUATORIAL_RADIUS + 300e3,
            0.001,
            98.0 |> deg2rad,
            ltdn_to_raan(10.5, jd₀),
            90.0 |> deg2rad,
            0.0
        )

        # Run a short analysis compiling the whole pipeline: the solver stack, the
        # right-hand side, and the output assembly.
        decay_analysis(
            orb;
            satellite_mass      = 100.0,
            satellite_mean_area = 1.0,
            gravity_model       = gm,
            space_indices       = (f107 = 140.0, f107_avg = 140.0, ap = 9.0),
            tf                  = 86400.0
        )

        # Compile the Jacchia 1977 and Jacchia-Roberts 1971 pipelines selected by the
        # macros `@decay_analysis__jacchia77` and `@decay_analysis__jr1971`. The space
        # indices provided by the macros are overridden with a constant named tuple to
        # keep the workload network-free.
        decay_analysis(
            orb;
            satellite_mass      = 100.0,
            satellite_mean_area = 1.0,
            gravity_model       = gm,
            tf                  = 86400.0,
            @decay_analysis__jacchia77,
            space_indices       = (f107 = 140.0, f107_avg = 140.0, kp = 3.0)
        )

        decay_analysis(
            orb;
            satellite_mass      = 100.0,
            satellite_mean_area = 1.0,
            gravity_model       = gm,
            tf                  = 86400.0,
            @decay_analysis__jr1971,
            space_indices       = (f107 = 140.0, f107_avg = 140.0, kp = 3.0)
        )

        # Compile the progress interface rendering.
        buf      = IOBuffer()
        progress = DecayProgress(buf, 86400.0, 300e3, 120e3; ansi = true)

        _start_decay_progress!(progress, 310e3)
        _update_decay_progress!(progress, 43200.0, 200e3, 210e3)
        _finish_decay_progress!(progress, 86400.0, 120e3, 130e3, true)
        _cleanup_decay_progress!(progress)

        rm(path; force = true)
    end
end
