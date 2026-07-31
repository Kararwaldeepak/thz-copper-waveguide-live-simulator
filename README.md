# THz Circular Waveguide — Live Browser Simulator

An interactive, browser-only simulator for coupling a focused Gaussian THz beam
into an air-filled circular copper waveguide.

The default experiment matches:

- copper tube length: **40 cm**
- waveguide inner diameter: **2 mm**
- lens focal length: **100 mm**
- lens clear diameter: **2 inches (50.8 mm)**
- frequency range: **0.10–2.00 THz**

Open `index.html` and change any control. The transverse field, cutoff,
Gaussian waist, modal overlap, conductor loss, transmitted power, group delay,
and spectral response update immediately. No Python, server, build step, or
installed software is required.

## Publish with GitHub Pages

1. Upload these files to the root of your GitHub repository:
   `index.html`, `styles.css`, `app.js`, `.nojekyll`, `README.md`, and
   `LICENSE`.
2. In the repository, open **Settings → Pages**.
3. Under **Build and deployment**, select **Deploy from a branch**.
4. Select the `main` branch and the `/ (root)` folder, then click **Save**.

For the repository `Kararwaldeepak/thz-circular-waveguide-modes`, the published
address will be:

<https://kararwaldeepak.github.io/thz-circular-waveguide-modes/>

GitHub may take a few minutes to publish the first deployment.

## What the browser calculates

- TE and TM circular-waveguide eigenfields from Bessel functions
- cutoff frequency: `fc = cχ / (2πa)`
- propagation constant: `β = √(k0² − kc²)`
- diffraction-limited Gaussian waist: `w0 = λf / (πwL)`
- centered x-polarized Gaussian-to-mode overlap
- first-order copper conductor loss using surface resistance
- output power after the selected copper-tube length
- relative group delay compared with free space

Available modes: TE11, TM01, TE21, TE01, TM11, and TE12.

## Scope

This is a fast analytical mode simulator, not a full 3D FDTD/FEM solver. It
uses an ideal circular air-filled guide with copper surface-resistance loss.
Flanges, wall roughness, beam offset or tilt, lens aberration, material
dispersion, and mode conversion along the tube are not included. Results are
suited to design exploration and should be checked against measurement or a
validated full-wave solver before fabrication.

## License

MIT License — see `LICENSE`.
