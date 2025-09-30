import io
import numpy as np
import streamlit as st

# Astropy pieces for units, coords, and cosmology-free geometry
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, AltAz, FK5, Galactic, GeocentricTrueEcliptic
from astropy.time import Time

st.set_page_config(page_title="Astro Basics Toolkit", page_icon="🛰️", layout="wide")

st.title("🛰️ Astro Basics Toolkit")
st.caption("RA/Dec • Frame transforms • Separation • Angular ↔ Physical • Pixel scale & FOV")
st.markdown("""
**Creator:** Anumanchi Agastya Sai Ram Likhit  
**GitHub:** [https://github.com/astropi-b/](https://github.com/astropi-b/)
""")

with st.sidebar:
    st.header("Quick Notes")
    st.markdown("""
- **Angles**: 1° = 60′ = 3600″  
- **Small-angle** (no cosmology): *size* = *θ* (radians) × *distance*  
- **Pixel scale**: arcsec/pixel = FOV(arcsec) / Npix  
- Use **ICRS** (J2000) for RA/Dec unless your data say otherwise.
""")

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "RA/Dec Converter", "Frame Transforms & Separation",
    "Angular ↔ Physical Size", "Pixel Scale & FOV", "Unit Snippets"
])

# =========================================================
# TAB 1 — RA/Dec Converter
# =========================================================
with tab1:
    st.subheader("RA/Dec Converter (sexagesimal ↔ decimal)")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Input (any format Astropy understands):**")
        ra_in = st.text_input("RA (e.g., '12:34:56.7' or '188.736°')", "12:34:56.7")
        dec_in = st.text_input("Dec (e.g., '-01:23:45.6' or '-1.395°')", "-01:23:45.6")
        frame_pick = st.selectbox("Frame", ["icrs", "fk5 J2000"], index=0)
        if frame_pick.startswith("fk5"):
            frame_obj = FK5(equinox=Time("J2000"))
        else:
            frame_obj = "icrs"

        if st.button("Convert"):
            try:
                c = SkyCoord(ra_in, dec_in, frame=frame_obj, unit=(u.hourangle, u.deg), obstime=Time("J2000"))
            except Exception:
                # Try degree-degree
                try:
                    c = SkyCoord(ra_in, dec_in, frame=frame_obj, unit=(u.deg, u.deg), obstime=Time("J2000"))
                except Exception as e:
                    st.error(f"Could not parse inputs: {e}")
                    c = None

            if c is not None:
                with c2:
                    st.markdown("**Outputs**")
                    st.code(
                        f"RA (decimal deg) : {c.ra.deg:.8f}\n"
                        f"Dec (decimal deg): {c.dec.deg:.8f}\n\n"
                        f"RA (HMS)         : {c.ra.to_string(unit=u.hour, sep=':')}\n"
                        f"Dec (DMS)        : {c.dec.to_string(unit=u.deg,  sep=':')}",
                        language="text",
                    )

# =========================================================
# TAB 2 — Frame Transforms & Separation
# =========================================================
with tab2:
    st.subheader("Coordinate Frame Transforms (ICRS ↔ Galactic ↔ Ecliptic)")
    colA, colB = st.columns(2)
    with colA:
        ra_t = st.text_input("RA (ICRS) for transform", "10:00:00")
        dec_t = st.text_input("Dec (ICRS) for transform", "+20:00:00")
        do_transform = st.button("Transform Frames")

    with colB:
        if do_transform:
            try:
                ic = SkyCoord(ra_t, dec_t, unit=(u.hourangle, u.deg), frame="icrs")
            except Exception:
                try:
                    ic = SkyCoord(ra_t, dec_t, unit=(u.deg, u.deg), frame="icrs")
                except Exception as e:
                    st.error(f"Parse error: {e}")
                    ic = None

        else:
            ic = None

        if ic is not None:
            gal = ic.transform_to(Galactic())
            ecl = ic.transform_to(GeocentricTrueEcliptic(obstime=Time("J2000")))

            st.markdown("**ICRS (J2000)**")
            st.code(f"RA={ic.ra.deg:.6f}°, Dec={ic.dec.deg:.6f}°", language="text")

            st.markdown("**Galactic**")
            st.code(f"ℓ={gal.l.deg:.6f}°, b={gal.b.deg:.6f}°", language="text")

            st.markdown("**Ecliptic (true, geocentric)**")
            st.code(f"λ={ecl.lon.deg:.6f}°, β={ecl.lat.deg:.6f}°", language="text")

    st.markdown("---")
    st.subheader("Great-Circle Separation")
    c1, c2 = st.columns(2)
    with c1:
        ra1 = st.text_input("RA1", "12:00:00")
        dec1 = st.text_input("Dec1", "+15:00:00")
        ra2 = st.text_input("RA2", "12:30:00")
        dec2 = st.text_input("Dec2", "+16:00:00")
        if st.button("Compute Separation"):
            try:
                a = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
            except Exception:
                a = SkyCoord(ra1, dec1, unit=(u.deg, u.deg))
            try:
                b = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg))
            except Exception:
                b = SkyCoord(ra2, dec2, unit=(u.deg, u.deg))
            sep = a.separation(b)
            with c2:
                st.code(
                    f"Separation: {sep.deg:.6f}°  = {sep.arcmin:.3f}′  = {sep.arcsec:.3f}″",
                    language="text",
                )

# =========================================================
# TAB 3 — Angular ↔ Physical Size (Small-angle)
# =========================================================
with tab3:
    st.subheader("Angular ↔ Physical Size (no cosmology, small-angle)")
    col1, col2 = st.columns(2)

    with col1:
        mode = st.radio("Mode", ["Given distance: angle → physical", "Given distance: physical → angle"], index=0)
        distance_value = st.number_input("Distance", value=140.0, min_value=0.0, help="e.g., 140 pc to a nearby molecular cloud")
        distance_unit = st.selectbox("Distance unit", ["pc", "kpc", "Mpc", "ly"], index=0)

        if mode.startswith("Given distance: angle"):
            theta_input = st.number_input("Angular size", value=60.0, help="Arcseconds or arcminutes or degrees")
            theta_unit = st.selectbox("Angle unit", ["arcsec", "arcmin", "deg"], index=0)
        else:
            size_input = st.number_input("Physical size", value=0.1, help="e.g., 0.1 pc")
            size_unit = st.selectbox("Physical unit", ["m", "km", "AU", "pc", "kpc", "ly"], index=3)

        if st.button("Compute (small-angle)"):
            # Units
            D = (distance_value * getattr(u, distance_unit))
            if mode.startswith("Given distance: angle"):
                theta = (theta_input * getattr(u, theta_unit)).to(u.rad)
                phys = (D * theta).to(u.m)
                # Present in common astro units
                st.success("Physical size")
                st.code(
                    f"{phys.to(u.AU):.6g}  (AU)\n"
                    f"{phys.to(u.pc):.6g}  (pc)\n"
                    f"{phys.to(u.lyr):.6g}  (ly)",
                    language="text",
                )
            else:
                L = (size_input * getattr(u, size_unit)).to(u.m)
                theta = (L / D).to(u.rad)
                st.success("Angular size")
                st.code(
                    f"{theta.to(u.deg):.6g}  (deg)\n"
                    f"{theta.to(u.arcmin):.6g}  (arcmin)\n"
                    f"{theta.to(u.arcsec):.6g}  (arcsec)",
                    language="text",
                )

    with col2:
        st.info("Tip: For cosmological distances you’ll want **angular-diameter distance** (needs a cosmology). "
                "This tool intentionally uses **pure geometry** for nearby objects (e.g., ISM filaments).")

# =========================================================
# TAB 4 — Pixel Scale & FOV
# =========================================================
with tab4:
    st.subheader("Pixel Scale ⇄ FOV ⇄ Image Size")
    st.markdown("Supply **any two**, compute the third. All angles on sky.")

    col1, col2, col3 = st.columns(3)
    # Inputs
    with col1:
        fov_arcmin = st.number_input("FOV (arcmin)", value=10.0, min_value=0.0)
        fov_axis = st.selectbox("FOV axis", ["width", "height"], index=0)
    with col2:
        n_pixels = st.number_input("Pixels along axis (N)", value=2048, min_value=1, step=1)
    with col3:
        pix_scale_arcsec = st.number_input("Pixel scale (arcsec/pixel)", value=0.0, min_value=0.0)

    if st.button("Solve"):
        fov = (fov_arcmin * u.arcmin).to(u.arcsec).value  # arcsec
        N = float(n_pixels)
        s = pix_scale_arcsec  # arcsec/pix

        # Determine which is zero or missing -> compute it
        eps = 1e-12
        known = sum([fov > eps, N > eps, s > eps])
        if known < 2:
            st.error("Provide at least two values.")
        else:
            if s <= eps:  # compute pixel scale
                s = fov / N
            elif fov <= eps:  # compute FOV
                fov = s * N
            elif N <= eps:
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
        scale_mode = st.radio("Mode", ["Pixel scale → linear per pixel", "Linear per pixel → pixel scale"], index=0)
        D_val = st.number_input("Distance", value=140.0, min_value=0.0)
        D_unit = st.selectbox("Distance unit", ["pc", "kpc", "Mpc", "ly"], index=0)

        if scale_mode.startswith("Pixel scale"):
            s_arcsec = st.number_input("Pixel scale (arcsec/pixel)", value=1.0, min_value=0.0)
        else:
            L_m = st.number_input("Linear scale per pixel", value=50.0, min_value=0.0, help="e.g., AU per pixel if you pick AU below")
            L_unit = st.selectbox("Linear unit per pixel", ["m", "km", "AU", "pc", "ly"], index=2)

    if st.button("Compute linear/angle scale"):
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
            st.code(
                f"{theta:.6g} arcsec/pixel",
                language="text",
            )

# =========================================================
# TAB 5 — Unit Snippets / Cheats
# =========================================================
with tab5:
    st.subheader("Quick Unit Conversions (Astropy Units)")
    col1, col2, col3 = st.columns(3)
    with col1:
        ang_val = st.number_input("Angle value", value=1.0)
        ang_unit = st.selectbox("Angle unit", ["deg", "arcmin", "arcsec", "rad"], index=0)
        if st.button("Convert angle"):
            val = (ang_val * getattr(u, ang_unit))
            st.code(
                f"{val.to(u.deg):.8g} deg\n"
                f"{val.to(u.arcmin):.8g} arcmin\n"
                f"{val.to(u.arcsec):.8g} arcsec\n"
                f"{val.to(u.rad):.8g} rad",
                language="text",
            )
    with col2:
        len_val = st.number_input("Length value", value=1.0, key="lenv")
        len_unit = st.selectbox("Length unit", ["m", "km", "AU", "pc", "kpc", "Mpc", "ly"], index=3)
        if st.button("Convert length"):
            val = (len_val * getattr(u, len_unit))
            st.code(
                f"{val.to(u.m):.8g} m\n"
                f"{val.to(u.km):.8g} km\n"
                f"{val.to(u.AU):.8g} AU\n"
                f"{val.to(u.pc):.8g} pc\n"
                f"{val.to(u.kpc):.8g} kpc\n"
                f"{val.to(u.Mpc):.8g} Mpc\n"
                f"{val.to(u.lyr):.8g} ly",
                language="text",
            )
    with col3:
        st.info("This page intentionally avoids cosmology (z→DA) to keep it lightweight. "
                "We can add Planck18 later for angular-diameter distance if you want.")
