with tab4:
    st.subheader("Pixel Scale ⇄ FOV ⇄ Image Size")
    st.markdown("Supply **any two**, compute the third. All angles on sky.")

    col1, col2, col3 = st.columns(3)
    with col1:
        fov_arcmin = st.number_input("FOV (arcmin)", value=10.0, min_value=0.0, key="t4_fov_arcmin")
        fov_axis = st.selectbox("FOV axis", ["width", "height"], index=0, key="t4_fov_axis")
    with col2:
        n_pixels = st.number_input("Pixels along axis (N)", value=2048, min_value=1, step=1, key="t4_npix")
    with col3:
        pix_scale_arcsec = st.number_input("Pixel scale (arcsec/pixel)", value=0.0, min_value=0.0, key="t4_pixscale")

    if st.button("Solve", key="t4_solve"):
        import astropy.units as u
        fov = (fov_arcmin * u.arcmin).to(u.arcsec).value  # arcsec
        N = float(n_pixels)
        s = pix_scale_arcsec  # arcsec/pix

        eps = 1e-12
        known = sum([fov > eps, N > eps, s > eps])
        if known < 2:
            st.error("Provide at least two values.")
        else:
            if s <= eps:       # compute pixel scale
                s = fov / N
            elif fov <= eps:   # compute FOV
                fov = s * N
            elif N <= eps:     # compute N
                N = fov / s

            st.success("Results")
            st.code(
                f"FOV:  {fov/60:.6f} arcmin  = {fov/3600:.6f} deg\n"
                f"Npix: {int(round(N))}\n"
                f"Scale:{s:.6f} arcsec/pixel",
                language="text",
            )

    st.markdown("---")
    st.subheader("Angular ⇄ Linear scale at distance")
    colA, colB = st.columns(2)
    with colA:
        scale_mode = st.radio(
            "Mode",
            ["Pixel scale → linear per pixel", "Linear per pixel → pixel scale"],
            index=0, key="t4_scale_mode",
        )
        D_val = st.number_input("Distance", value=140.0, min_value=0.0, key="t4_D_val")
        D_unit = st.selectbox("Distance unit", ["pc", "kpc", "Mpc", "ly"], index=0, key="t4_D_unit")

        if scale_mode.startswith("Pixel scale"):
            s_arcsec = st.number_input("Pixel scale (arcsec/pixel)", value=1.0, min_value=0.0, key="t4_s_arcsec")
        else:
            L_m = st.number_input(
                "Linear scale per pixel", value=50.0, min_value=0.0,
                help="e.g., AU per pixel if you pick AU below",
                key="t4_L_m",
            )
            L_unit = st.selectbox("Linear unit per pixel", ["m", "km", "AU", "pc", "ly"], index=2, key="t4_L_unit")

    if st.button("Compute linear/angle scale", key="t4_compute_scale"):
        import astropy.units as u
        D = (D_val * getattr(u, D_unit))
        if scale_mode.startswith("Pixel scale"):
            theta = (s_arcsec * u.arcsec).to(u.rad)
            L = (D * theta).to(u.m)
            st.code(
                f"{L.to(u.AU):.6g} per pixel (AU)\n"
                f"{L.to(u.km):.6g} per pixel (km)\n"
                f"{L.to(u.m):.6g} per pixel (m)",
                language="text",
            )
        else:
            L = (L_m * getattr(u, L_unit)).to(u.m)
            theta = (L / D).to(u.arcsec)
            st.code(f"{theta:.6g} arcsec/pixel", language="text")
