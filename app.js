(() => {
  "use strict";

  const C0 = 299792458;
  const MU0 = 4 * Math.PI * 1e-7;
  const EPS0 = 1 / (MU0 * C0 * C0);
  const COPPER_SIGMA = 5.8e7;
  const DB_PER_NEPER = 8.685889638;
  const TWO_PI = 2 * Math.PI;

  const MODES = {
    TE11: { family: "TE", m: 1, n: 1, root: 1.841183781, orientation: "sin" },
    TM01: { family: "TM", m: 0, n: 1, root: 2.404825558, orientation: "cos" },
    TE21: { family: "TE", m: 2, n: 1, root: 3.054236928, orientation: "cos" },
    TE01: { family: "TE", m: 0, n: 1, root: 3.831705970, orientation: "cos" },
    TM11: { family: "TM", m: 1, n: 1, root: 3.831705970, orientation: "cos" },
    TE12: { family: "TE", m: 1, n: 2, root: 5.331442774, orientation: "sin" },
  };

  const DEFAULTS = {
    frequencyTHz: 1,
    diameterMM: 2,
    lengthCM: 40,
    focalMM: 100,
    lensDiameterMM: 50.8,
    fill: 1,
    mode: "TE11",
  };

  const CONTROL_MAP = {
    frequency: { key: "frequencyTHz", range: "frequencyRange", number: "frequencyNumber" },
    diameter: { key: "diameterMM", range: "diameterRange", number: "diameterNumber" },
    length: { key: "lengthCM", range: "lengthRange", number: "lengthNumber" },
    focal: { key: "focalMM", range: "focalRange", number: "focalNumber" },
    lensDiameter: {
      key: "lensDiameterMM",
      range: "lensDiameterRange",
      number: "lensDiameterNumber",
    },
    fill: { key: "fill", range: "fillRange", number: "fillNumber" },
  };

  const state = { ...DEFAULTS };
  const spatialCaches = new Map();
  const element = {};
  let renderQueued = false;
  let chartState = null;
  let chartDragging = false;

  function factorial(n) {
    let value = 1;
    for (let i = 2; i <= n; i += 1) value *= i;
    return value;
  }

  function besselJ(n, x) {
    let term = Math.pow(x / 2, n) / factorial(n);
    let sum = term;
    for (let k = 1; k < 60; k += 1) {
      term *= -(x * x * 0.25) / (k * (k + n));
      sum += term;
      if (Math.abs(term) < 1e-14 * Math.max(1, Math.abs(sum))) break;
    }
    return sum;
  }

  function besselPrime(n, x) {
    if (n === 0) return -besselJ(1, x);
    return 0.5 * (besselJ(n - 1, x) - besselJ(n + 1, x));
  }

  function angularTerms(mode, phi) {
    if (mode.m === 0) return { angular: 1, derivative: 0 };
    const angle = mode.m * phi;
    if (mode.orientation === "sin") {
      return {
        angular: Math.sin(angle),
        derivative: mode.m * Math.cos(angle),
      };
    }
    return {
      angular: Math.cos(angle),
      derivative: -mode.m * Math.sin(angle),
    };
  }

  function scalarAndGradient(x, y, mode) {
    const r = Math.hypot(x, y);
    if (r > 1.000001) return null;

    if (r < 1e-9) {
      if (mode.m === 0) {
        return { psi: 1, gx: 0, gy: 0 };
      }
      if (mode.m === 1) {
        const centerSlope = mode.root * 0.5;
        return mode.orientation === "sin"
          ? { psi: 0, gx: 0, gy: centerSlope }
          : { psi: 0, gx: centerSlope, gy: 0 };
      }
      return { psi: 0, gx: 0, gy: 0 };
    }

    const phi = Math.atan2(y, x);
    const cosPhi = x / r;
    const sinPhi = y / r;
    const { angular, derivative } = angularTerms(mode, phi);
    const radial = besselJ(mode.m, mode.root * r);
    const radialDerivative = mode.root * besselPrime(mode.m, mode.root * r) * angular;
    const angularDerivative = (radial * derivative) / r;

    return {
      psi: radial * angular,
      gx: radialDerivative * cosPhi - angularDerivative * sinPhi,
      gy: radialDerivative * sinPhi + angularDerivative * cosPhi,
    };
  }

  function getSpatialCache(modeKey) {
    if (spatialCaches.has(modeKey)) return spatialCaches.get(modeKey);

    const mode = MODES[modeKey];
    const resolution = 121;
    const step = 2 / (resolution - 1);
    const points = [];
    let gradientIntegral = 0;

    for (let row = 0; row < resolution; row += 1) {
      const y = 1 - row * step;
      for (let col = 0; col < resolution; col += 1) {
        const x = -1 + col * step;
        const sample = scalarAndGradient(x, y, mode);
        if (!sample) continue;

        const baseEx = mode.family === "TE" ? -sample.gy : -sample.gx;
        const baseEy = mode.family === "TE" ? sample.gx : -sample.gy;
        points.push({ row, col, x, y, ...sample, baseEx, baseEy });
        gradientIntegral += (sample.gx * sample.gx + sample.gy * sample.gy) * step * step;
      }
    }

    const cache = { resolution, step, points, gradientIntegral };
    spatialCaches.set(modeKey, cache);
    return cache;
  }

  function overlapEfficiency(modeKey, radius, waist) {
    const cache = getSpatialCache(modeKey);
    let overlap = 0;
    let gaussianNorm = 0;
    let modeNorm = 0;
    const areaWeight = cache.step * cache.step;

    for (const point of cache.points) {
      const physicalRadius = Math.hypot(point.x, point.y) * radius;
      const gaussian = Math.exp(-Math.pow(physicalRadius / waist, 2));
      overlap += gaussian * point.baseEx * areaWeight;
      gaussianNorm += gaussian * gaussian * areaWeight;
      modeNorm +=
        (point.baseEx * point.baseEx + point.baseEy * point.baseEy) * areaWeight;
    }

    const denominator = gaussianNorm * modeNorm;
    return denominator > 0 ? Math.min(1, Math.max(0, (overlap * overlap) / denominator)) : 0;
  }

  function conductorLoss(modeKey, frequencyHz, radius, beta) {
    if (!(beta > 0)) return { alpha: Number.NaN, dbPerM: Number.NaN };

    const mode = MODES[modeKey];
    const cache = getSpatialCache(modeKey);
    const omega = TWO_PI * frequencyHz;
    const boundarySamples = 720;
    const dPhi = TWO_PI / boundarySamples;
    let boundaryMagneticIntegral = 0;

    if (mode.family === "TE") {
      const eFactor = (omega * MU0 * radius) / (mode.root * mode.root);
      const hFactor = (beta * radius) / (mode.root * mode.root);
      const guidedPower =
        0.5 * eFactor * hFactor * radius * radius * cache.gradientIntegral;

      for (let i = 0; i < boundarySamples; i += 1) {
        const phi = (i + 0.5) * dPhi;
        const { angular, derivative } = angularTerms(mode, phi);
        const psiBoundary = besselJ(mode.m, mode.root) * angular;
        const tangentialGradient = besselJ(mode.m, mode.root) * derivative;
        const hPhi = hFactor * tangentialGradient;
        boundaryMagneticIntegral += hPhi * hPhi + psiBoundary * psiBoundary;
      }

      boundaryMagneticIntegral *= radius * dPhi;
      const surfaceResistance = Math.sqrt(
        (Math.PI * frequencyHz * MU0) / COPPER_SIGMA,
      );
      const wallLoss = 0.5 * surfaceResistance * boundaryMagneticIntegral;
      const alpha = wallLoss / (2 * guidedPower);
      return { alpha, dbPerM: DB_PER_NEPER * alpha };
    }

    const eFactor = (beta * radius) / (mode.root * mode.root);
    const hFactor = (omega * EPS0 * radius) / (mode.root * mode.root);
    const guidedPower =
      0.5 * eFactor * hFactor * radius * radius * cache.gradientIntegral;

    for (let i = 0; i < boundarySamples; i += 1) {
      const phi = (i + 0.5) * dPhi;
      const { angular } = angularTerms(mode, phi);
      const radialGradient = mode.root * besselPrime(mode.m, mode.root) * angular;
      const hPhi = hFactor * radialGradient;
      boundaryMagneticIntegral += hPhi * hPhi;
    }

    boundaryMagneticIntegral *= radius * dPhi;
    const surfaceResistance = Math.sqrt(
      (Math.PI * frequencyHz * MU0) / COPPER_SIGMA,
    );
    const wallLoss = 0.5 * surfaceResistance * boundaryMagneticIntegral;
    const alpha = wallLoss / (2 * guidedPower);
    return { alpha, dbPerM: DB_PER_NEPER * alpha };
  }

  function solve(parameters = state) {
    const mode = MODES[parameters.mode];
    const frequencyHz = parameters.frequencyTHz * 1e12;
    const radius = parameters.diameterMM * 0.5e-3;
    const length = parameters.lengthCM * 1e-2;
    const focalLength = parameters.focalMM * 1e-3;
    const lensDiameter = parameters.lensDiameterMM * 1e-3;
    const illuminatedLensRadius = parameters.fill * lensDiameter * 0.5;
    const wavelength = C0 / frequencyHz;
    const waist = (wavelength * focalLength) / (Math.PI * illuminatedLensRadius);
    const numericalAperture = Math.sin(Math.atan((lensDiameter * 0.5) / focalLength));
    const fNumber = focalLength / lensDiameter;
    const apertureThroughput = 1 - Math.exp(-2 * Math.pow(radius / waist, 2));
    const cutoffHz = (C0 * mode.root) / (TWO_PI * radius);
    const k0 = TWO_PI * frequencyHz / C0;
    const kc = mode.root / radius;
    const propagating = frequencyHz > cutoffHz;
    const beta = propagating ? Math.sqrt(k0 * k0 - kc * kc) : Number.NaN;

    const normalizedOverlap = overlapEfficiency(parameters.mode, radius, waist);
    const coupling = propagating
      ? Math.min(1, Math.max(0, apertureThroughput * normalizedOverlap))
      : 0;

    let alpha = Number.NaN;
    let lossDbPerM = Number.NaN;
    let totalLossDb = Number.NaN;
    let outputPower = 0;
    let relativeGroupDelayPs = Number.NaN;

    if (propagating) {
      const loss = conductorLoss(parameters.mode, frequencyHz, radius, beta);
      alpha = loss.alpha;
      lossDbPerM = loss.dbPerM;
      totalLossDb = lossDbPerM * length;
      outputPower = coupling * Math.exp(-2 * alpha * length);
      const groupVelocity = C0 * Math.sqrt(1 - Math.pow(cutoffHz / frequencyHz, 2));
      relativeGroupDelayPs = length * (1 / groupVelocity - 1 / C0) * 1e12;
    }

    return {
      mode,
      frequencyHz,
      radius,
      length,
      wavelength,
      waist,
      numericalAperture,
      fNumber,
      apertureThroughput,
      cutoffHz,
      k0,
      kc,
      propagating,
      beta,
      normalizedOverlap,
      coupling,
      alpha,
      lossDbPerM,
      totalLossDb,
      outputPower,
      relativeGroupDelayPs,
    };
  }

  function clamp(value, min, max) {
    return Math.min(max, Math.max(min, value));
  }

  function formatModeHtml(modeKey) {
    const mode = MODES[modeKey];
    return `${mode.family}<sub>${mode.m}${mode.n}</sub>`;
  }

  function setText(id, value) {
    element[id].textContent = value;
  }

  function formatFinite(value, digits) {
    return Number.isFinite(value) ? value.toFixed(digits) : "—";
  }

  function updateResults(result, elapsedMs) {
    element.modeLabel.innerHTML = formatModeHtml(state.mode);
    setText("frequencyLabel", `${state.frequencyTHz.toFixed(2)} THz`);
    setText("cutoffValue", (result.cutoffHz / 1e9).toFixed(2));
    element.cutoffFormula.innerHTML = `f<sub>c</sub>(${formatModeHtml(state.mode)})`;
    setText("waistValue", (result.waist * 1e3).toFixed(3));
    setText("naValue", result.numericalAperture.toFixed(3));
    setText("fNumberValue", `f/${result.fNumber.toFixed(2)}`);
    setText("couplingValue", result.propagating ? (result.coupling * 100).toFixed(1) : "—");
    setText(
      "transmissionValue",
      result.propagating ? (result.outputPower * 100).toFixed(1) : "—",
    );
    setText("lossValue", formatFinite(result.totalLossDb, 2));
    setText(
      "lossPerMeterValue",
      Number.isFinite(result.lossDbPerM)
        ? `${result.lossDbPerM.toFixed(2)} dB m⁻¹`
        : "— dB m⁻¹",
    );
    setText("delayValue", formatFinite(result.relativeGroupDelayPs, 2));
    setText("rootValue", result.mode.root.toFixed(5));
    setText(
      "betaValue",
      result.propagating ? `${(result.beta / 1e3).toFixed(2)} mm⁻¹` : "imaginary",
    );
    setText("apertureValue", `${(result.apertureThroughput * 100).toFixed(1)}%`);
    setText("solverTime", `${elapsedMs.toFixed(1)} ms`);

    element.fieldStatus.textContent = result.propagating ? "Propagating" : "Evanescent";
    element.fieldStatus.classList.toggle("evanescent", !result.propagating);
    element.statusCard.classList.toggle("evanescent", !result.propagating);
    setText("propagationStatus", result.propagating ? "Propagating" : "Below cutoff");
    setText(
      "propagationDetail",
      result.propagating
        ? "β is real above cutoff"
        : `Raise f above ${(result.cutoffHz / 1e12).toFixed(3)} THz`,
    );
  }

  function resizeCanvas(canvas) {
    const rect = canvas.getBoundingClientRect();
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const width = Math.max(1, Math.round(rect.width * dpr));
    const height = Math.max(1, Math.round(rect.height * dpr));
    if (canvas.width !== width || canvas.height !== height) {
      canvas.width = width;
      canvas.height = height;
    }
    const context = canvas.getContext("2d");
    context.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { context, width: rect.width, height: rect.height, dpr };
  }

  function hexToRgb(hex) {
    const value = Number.parseInt(hex.slice(1), 16);
    return {
      r: (value >> 16) & 255,
      g: (value >> 8) & 255,
      b: value & 255,
    };
  }

  function interpolateColor(first, second, t) {
    const a = hexToRgb(first);
    const b = hexToRgb(second);
    const mix = (left, right) => Math.round(left + (right - left) * t);
    return `rgb(${mix(a.r, b.r)}, ${mix(a.g, b.g)}, ${mix(a.b, b.b)})`;
  }

  function fieldColor(value, dimmed) {
    const stops = [
      [-1, "#071e57"],
      [-0.48, "#1454ba"],
      [0, "#07131f"],
      [0.48, "#1dd6da"],
      [0.76, "#f6b73c"],
      [1, "#ff7043"],
    ];
    const v = clamp(value, -1, 1);
    let color = stops[stops.length - 1][1];
    for (let i = 0; i < stops.length - 1; i += 1) {
      const [leftValue, leftColor] = stops[i];
      const [rightValue, rightColor] = stops[i + 1];
      if (v >= leftValue && v <= rightValue) {
        color = interpolateColor(
          leftColor,
          rightColor,
          (v - leftValue) / (rightValue - leftValue),
        );
        break;
      }
    }
    if (!dimmed) return color;
    const rgb = color.match(/\d+/g).map(Number);
    return `rgba(${rgb[0]}, ${rgb[1]}, ${rgb[2]}, 0.38)`;
  }

  function drawArrow(context, x, y, dx, dy, scale, alpha) {
    const magnitude = Math.hypot(dx, dy);
    if (magnitude < 1e-6) return;
    const ux = dx / magnitude;
    const uy = -dy / magnitude;
    const length = scale * (0.32 + 0.68 * magnitude);
    const half = length * 0.5;
    const startX = x - ux * half;
    const startY = y - uy * half;
    const endX = x + ux * half;
    const endY = y + uy * half;
    const head = Math.max(3, length * 0.22);

    context.beginPath();
    context.moveTo(startX, startY);
    context.lineTo(endX, endY);
    context.moveTo(endX, endY);
    context.lineTo(endX - ux * head + uy * head * 0.55, endY - uy * head - ux * head * 0.55);
    context.moveTo(endX, endY);
    context.lineTo(endX - ux * head - uy * head * 0.55, endY - uy * head + ux * head * 0.55);
    context.strokeStyle = `rgba(225, 250, 255, ${alpha})`;
    context.lineWidth = 0.9;
    context.stroke();
  }

  function drawField(result) {
    const { context, width, height } = resizeCanvas(element.fieldCanvas);
    const cache = getSpatialCache(state.mode);
    const mode = MODES[state.mode];
    const padding = 18;
    const diameter = Math.max(80, Math.min(width, height) - padding * 2);
    const left = (width - diameter) * 0.5;
    const top = (height - diameter) * 0.5;
    const cell = diameter / cache.resolution;
    const betaScale =
      mode.family === "TM" && result.propagating
        ? result.beta / result.kc
        : 1;

    context.clearRect(0, 0, width, height);
    context.fillStyle = "#050c14";
    context.fillRect(0, 0, width, height);

    let maxField = 0;
    for (const point of cache.points) {
      maxField = Math.max(maxField, Math.abs(point.baseEx * betaScale));
    }
    maxField = maxField || 1;

    for (const point of cache.points) {
      const value = (point.baseEx * betaScale) / maxField;
      context.fillStyle = fieldColor(value, !result.propagating);
      context.fillRect(
        left + point.col * cell,
        top + point.row * cell,
        cell + 0.75,
        cell + 0.75,
      );
    }

    context.save();
    context.beginPath();
    context.arc(width * 0.5, height * 0.5, diameter * 0.5, 0, TWO_PI);
    context.clip();
    context.setLineDash([5, 6]);
    context.strokeStyle = "rgba(160, 192, 208, 0.35)";
    context.lineWidth = 1;
    context.beginPath();
    context.moveTo(left, height * 0.5);
    context.lineTo(left + diameter, height * 0.5);
    context.moveTo(width * 0.5, top);
    context.lineTo(width * 0.5, top + diameter);
    context.stroke();
    context.setLineDash([]);

    const arrowStep = 15;
    for (let row = 7; row < cache.resolution; row += arrowStep) {
      for (let col = 7; col < cache.resolution; col += arrowStep) {
        const xUnit = -1 + col * cache.step;
        const yUnit = 1 - row * cache.step;
        const sample = scalarAndGradient(xUnit, yUnit, mode);
        if (!sample) continue;
        let ex = mode.family === "TE" ? -sample.gy : -sample.gx;
        let ey = mode.family === "TE" ? sample.gx : -sample.gy;
        ex *= betaScale;
        ey *= betaScale;
        const normalizedMagnitude = clamp(Math.hypot(ex, ey) / maxField, 0, 1);
        drawArrow(
          context,
          left + col * cell,
          top + row * cell,
          ex / maxField,
          ey / maxField,
          cell * arrowStep * 0.72,
          result.propagating ? 0.25 + 0.65 * normalizedMagnitude : 0.18,
        );
      }
    }
    context.restore();

    context.beginPath();
    context.arc(width * 0.5, height * 0.5, diameter * 0.5 + 1, 0, TWO_PI);
    context.strokeStyle = result.propagating ? "#f5a524" : "#ff5577";
    context.lineWidth = 2;
    context.shadowColor = result.propagating
      ? "rgba(245, 165, 36, 0.48)"
      : "rgba(255, 85, 119, 0.42)";
    context.shadowBlur = 9;
    context.stroke();
    context.shadowBlur = 0;
  }

  function frequencySweep() {
    const points = [];
    const count = 121;
    for (let i = 0; i < count; i += 1) {
      const frequencyTHz = 0.1 + (1.9 * i) / (count - 1);
      const result = solve({ ...state, frequencyTHz });
      points.push({
        frequencyTHz,
        outputPercent: result.outputPower * 100,
        lossDb: result.totalLossDb,
        propagating: result.propagating,
      });
    }
    return points;
  }

  function drawSpectrum(currentResult) {
    const { context, width, height } = resizeCanvas(element.spectrumCanvas);
    const margin = { left: 44, right: 42, top: 13, bottom: 26 };
    const plot = {
      left: margin.left,
      top: margin.top,
      width: Math.max(10, width - margin.left - margin.right),
      height: Math.max(10, height - margin.top - margin.bottom),
    };
    const points = frequencySweep();
    const finiteLosses = points.filter((point) => Number.isFinite(point.lossDb)).map((point) => point.lossDb);
    const rawLossMax = Math.max(1, ...finiteLosses);
    const lossMax = Math.ceil(rawLossMax / 2) * 2;
    const xForFrequency = (frequencyTHz) =>
      plot.left + ((frequencyTHz - 0.1) / 1.9) * plot.width;
    const yForOutput = (value) => plot.top + plot.height * (1 - clamp(value, 0, 100) / 100);
    const yForLoss = (value) => plot.top + plot.height * (clamp(value, 0, lossMax) / lossMax);

    context.clearRect(0, 0, width, height);
    context.fillStyle = "#07111c";
    context.fillRect(0, 0, width, height);
    context.font = '9px "SFMono-Regular", Consolas, monospace';
    context.textBaseline = "middle";

    for (let i = 0; i <= 4; i += 1) {
      const y = plot.top + (plot.height * i) / 4;
      context.beginPath();
      context.moveTo(plot.left, y);
      context.lineTo(plot.left + plot.width, y);
      context.strokeStyle = "rgba(74, 113, 139, 0.23)";
      context.lineWidth = 1;
      context.stroke();
      context.fillStyle = "#6f879a";
      context.textAlign = "right";
      context.fillText(`${100 - i * 25}`, plot.left - 6, y);
      context.textAlign = "left";
      context.fillText(`${((lossMax * i) / 4).toFixed(1)}`, plot.left + plot.width + 6, y);
    }

    for (let i = 0; i <= 4; i += 1) {
      const frequency = 0.1 + (1.9 * i) / 4;
      const x = xForFrequency(frequency);
      context.beginPath();
      context.moveTo(x, plot.top);
      context.lineTo(x, plot.top + plot.height);
      context.strokeStyle = "rgba(74, 113, 139, 0.16)";
      context.stroke();
      context.fillStyle = "#6f879a";
      context.textAlign = "center";
      context.fillText(frequency.toFixed(2), x, plot.top + plot.height + 13);
    }

    context.save();
    context.beginPath();
    context.rect(plot.left, plot.top, plot.width, plot.height);
    context.clip();

    const cutoffTHz = currentResult.cutoffHz / 1e12;
    if (cutoffTHz >= 0.1 && cutoffTHz <= 2) {
      const cutoffX = xForFrequency(cutoffTHz);
      context.beginPath();
      context.moveTo(cutoffX, plot.top);
      context.lineTo(cutoffX, plot.top + plot.height);
      context.setLineDash([4, 4]);
      context.strokeStyle = "rgba(184, 204, 215, 0.72)";
      context.lineWidth = 1;
      context.stroke();
      context.setLineDash([]);
    }

    context.beginPath();
    let outputStarted = false;
    for (const point of points) {
      if (!point.propagating) continue;
      const x = xForFrequency(point.frequencyTHz);
      const y = yForOutput(point.outputPercent);
      if (!outputStarted) {
        context.moveTo(x, y);
        outputStarted = true;
      } else {
        context.lineTo(x, y);
      }
    }
    context.strokeStyle = "#32d7e7";
    context.lineWidth = 2;
    context.shadowColor = "rgba(50, 215, 231, 0.35)";
    context.shadowBlur = 7;
    context.stroke();
    context.shadowBlur = 0;

    context.beginPath();
    let lossStarted = false;
    for (const point of points) {
      if (!Number.isFinite(point.lossDb)) continue;
      const x = xForFrequency(point.frequencyTHz);
      const y = yForLoss(point.lossDb);
      if (!lossStarted) {
        context.moveTo(x, y);
        lossStarted = true;
      } else {
        context.lineTo(x, y);
      }
    }
    context.strokeStyle = "#f5a524";
    context.lineWidth = 1.8;
    context.stroke();

    const currentX = xForFrequency(state.frequencyTHz);
    context.beginPath();
    context.moveTo(currentX, plot.top);
    context.lineTo(currentX, plot.top + plot.height);
    context.strokeStyle = "rgba(50, 215, 231, 0.55)";
    context.lineWidth = 1;
    context.stroke();
    if (currentResult.propagating) {
      context.beginPath();
      context.arc(currentX, yForOutput(currentResult.outputPower * 100), 3.5, 0, TWO_PI);
      context.fillStyle = "#e8fdff";
      context.fill();
      context.strokeStyle = "#32d7e7";
      context.lineWidth = 2;
      context.stroke();
    }
    context.restore();

    context.fillStyle = "#8aa0b4";
    context.textAlign = "left";
    context.fillText("output %", plot.left, 5);
    context.textAlign = "right";
    context.fillText("loss dB", plot.left + plot.width, 5);
    context.textAlign = "center";
    context.fillText("frequency (THz)", plot.left + plot.width * 0.5, height - 5);

    chartState = { points, plot, xForFrequency, lossMax };
  }

  function updateRangeFill(input) {
    const min = Number(input.min);
    const max = Number(input.max);
    const fraction = (Number(input.value) - min) / (max - min);
    input.style.setProperty("--fill", `${clamp(fraction, 0, 1) * 100}%`);
  }

  function syncControlsFromState() {
    for (const config of Object.values(CONTROL_MAP)) {
      const range = element[config.range];
      const number = element[config.number];
      range.value = state[config.key];
      number.value = state[config.key];
      updateRangeFill(range);
    }
    for (const button of element.modeButtons.querySelectorAll(".mode-button")) {
      button.classList.toggle("active", button.dataset.mode === state.mode);
    }
  }

  function scheduleRender() {
    if (renderQueued) return;
    renderQueued = true;
    requestAnimationFrame(() => {
      renderQueued = false;
      render();
    });
  }

  function render() {
    const started = performance.now();
    const result = solve();
    drawField(result);
    drawSpectrum(result);
    updateResults(result, performance.now() - started);
  }

  function setFrequencyFromChart(clientX) {
    if (!chartState) return;
    const rect = element.spectrumCanvas.getBoundingClientRect();
    const localX = clientX - rect.left;
    const fraction = clamp(
      (localX - chartState.plot.left) / chartState.plot.width,
      0,
      1,
    );
    state.frequencyTHz = Math.round((0.1 + 1.9 * fraction) * 100) / 100;
    element.frequencyRange.value = state.frequencyTHz;
    element.frequencyNumber.value = state.frequencyTHz.toFixed(2);
    updateRangeFill(element.frequencyRange);
    scheduleRender();
  }

  function updateChartTooltip(event) {
    if (!chartState || chartDragging) {
      element.chartTooltip.hidden = true;
      return;
    }
    const rect = element.spectrumCanvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    const { plot, points } = chartState;
    if (
      x < plot.left ||
      x > plot.left + plot.width ||
      y < plot.top ||
      y > plot.top + plot.height
    ) {
      element.chartTooltip.hidden = true;
      return;
    }
    const frequency = 0.1 + 1.9 * ((x - plot.left) / plot.width);
    const index = clamp(
      Math.round(((frequency - 0.1) / 1.9) * (points.length - 1)),
      0,
      points.length - 1,
    );
    const point = points[index];
    element.chartTooltip.innerHTML = point.propagating
      ? `<b>${point.frequencyTHz.toFixed(3)} THz</b><br>Output ${point.outputPercent.toFixed(2)}%<br>Loss ${point.lossDb.toFixed(3)} dB`
      : `<b>${point.frequencyTHz.toFixed(3)} THz</b><br>Below cutoff`;
    element.chartTooltip.hidden = false;
    const tooltipLeft = clamp(x + 12, 6, rect.width - 126);
    const tooltipTop = clamp(y - 24, 4, rect.height - 54);
    element.chartTooltip.style.left = `${tooltipLeft}px`;
    element.chartTooltip.style.top = `${tooltipTop}px`;
  }

  function csvValue(value, digits = 8) {
    return Number.isFinite(value) ? value.toFixed(digits) : "";
  }

  function exportCsv() {
    const rows = [
      [
        "frequency_thz",
        "mode",
        "cutoff_thz",
        "waist_mm",
        "aperture_throughput",
        "modal_coupling",
        "loss_db_per_m",
        "total_loss_db",
        "output_power_fraction",
        "relative_group_delay_ps",
      ],
    ];

    for (let i = 0; i <= 190; i += 1) {
      const frequencyTHz = 0.1 + i * 0.01;
      const result = solve({ ...state, frequencyTHz });
      rows.push([
        frequencyTHz.toFixed(2),
        state.mode,
        (result.cutoffHz / 1e12).toFixed(8),
        (result.waist * 1e3).toFixed(8),
        result.apertureThroughput.toFixed(8),
        result.coupling.toFixed(8),
        csvValue(result.lossDbPerM),
        csvValue(result.totalLossDb),
        result.outputPower.toFixed(8),
        csvValue(result.relativeGroupDelayPs),
      ]);
    }

    const blob = new Blob([rows.map((row) => row.join(",")).join("\n")], {
      type: "text/csv;charset=utf-8",
    });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `thz-waveguide-${state.mode.toLowerCase()}-spectrum.csv`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }

  function bindControls() {
    for (const config of Object.values(CONTROL_MAP)) {
      const range = element[config.range];
      const number = element[config.number];
      const min = Number(range.min);
      const max = Number(range.max);

      range.addEventListener("input", () => {
        state[config.key] = Number(range.value);
        number.value = range.value;
        updateRangeFill(range);
        scheduleRender();
      });

      const updateFromNumber = () => {
        if (number.value === "") return;
        const value = clamp(Number(number.value), min, max);
        if (!Number.isFinite(value)) return;
        state[config.key] = value;
        range.value = value;
        updateRangeFill(range);
        scheduleRender();
      };
      number.addEventListener("input", updateFromNumber);
      number.addEventListener("change", () => {
        updateFromNumber();
        number.value = state[config.key];
      });
    }

    element.modeButtons.addEventListener("click", (event) => {
      const button = event.target.closest(".mode-button");
      if (!button) return;
      state.mode = button.dataset.mode;
      syncControlsFromState();
      scheduleRender();
    });

    element.resetButton.addEventListener("click", () => {
      Object.assign(state, DEFAULTS);
      syncControlsFromState();
      scheduleRender();
    });
    element.exportButton.addEventListener("click", exportCsv);

    element.spectrumCanvas.addEventListener("pointerdown", (event) => {
      chartDragging = true;
      element.spectrumCanvas.setPointerCapture(event.pointerId);
      element.chartTooltip.hidden = true;
      setFrequencyFromChart(event.clientX);
    });
    element.spectrumCanvas.addEventListener("pointermove", (event) => {
      if (chartDragging) setFrequencyFromChart(event.clientX);
      else updateChartTooltip(event);
    });
    element.spectrumCanvas.addEventListener("pointerup", (event) => {
      chartDragging = false;
      if (element.spectrumCanvas.hasPointerCapture(event.pointerId)) {
        element.spectrumCanvas.releasePointerCapture(event.pointerId);
      }
    });
    element.spectrumCanvas.addEventListener("pointercancel", () => {
      chartDragging = false;
    });
    element.spectrumCanvas.addEventListener("pointerleave", () => {
      if (!chartDragging) element.chartTooltip.hidden = true;
    });

    const resizeObserver = new ResizeObserver(scheduleRender);
    resizeObserver.observe(element.fieldCanvas);
    resizeObserver.observe(element.spectrumCanvas);
  }

  function initialize() {
    const ids = [
      "frequencyRange",
      "frequencyNumber",
      "diameterRange",
      "diameterNumber",
      "lengthRange",
      "lengthNumber",
      "focalRange",
      "focalNumber",
      "lensDiameterRange",
      "lensDiameterNumber",
      "fillRange",
      "fillNumber",
      "modeButtons",
      "resetButton",
      "exportButton",
      "fieldCanvas",
      "spectrumCanvas",
      "chartTooltip",
      "modeLabel",
      "frequencyLabel",
      "fieldStatus",
      "cutoffValue",
      "cutoffFormula",
      "waistValue",
      "naValue",
      "fNumberValue",
      "couplingValue",
      "transmissionValue",
      "lossValue",
      "lossPerMeterValue",
      "delayValue",
      "statusCard",
      "propagationStatus",
      "propagationDetail",
      "rootValue",
      "betaValue",
      "apertureValue",
      "solverTime",
    ];
    for (const id of ids) element[id] = document.getElementById(id);
    syncControlsFromState();
    bindControls();
    scheduleRender();
  }

  initialize();
})();
