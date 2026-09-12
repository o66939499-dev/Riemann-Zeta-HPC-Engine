function riemannTheta(t) {
  const term1 = (t / 2) * Math.log(t / (2 * Math.PI));
  const term2 = t / 2;
  const term3 = Math.PI / 8;
  const term4 = 1 / (48 * t);
  return term1 - term2 - term3 + term4;
}

function Z(t) {
  const theta = riemannTheta(t);
  const terms = Math.floor(Math.sqrt(t / (2 * Math.PI))) + 2;
  let realSum = 0;
  let imagSum = 0;
  for (let k = 1; k <= terms; k++) {
    const sign = (k % 2 === 0) ? -1 : 1;
    const mag = sign / Math.sqrt(k);
    const arg = -t * Math.log(k);
    realSum += mag * Math.cos(arg);
    imagSum += mag * Math.sin(arg);
  }
  return realSum * Math.cos(theta) - imagSum * Math.sin(theta);
}

// Brent's method: inverse quadratic interpolation + secant + bisection.
function brentZero(a, b, tol = 1e-10, maxIter = 100) {
  let fa = Z(a);
  let fb = Z(b);
  if (fa * fb >= 0) return (a + b) / 2;

  if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }

  let c = a, fc = fa, d = b - a, mflag = true;

  for (let i = 0; i < maxIter; i++) {
    if (fb === 0 || Math.abs(b - a) < tol) break;

    let s;
    if (fa !== fc && fb !== fc) {
      s = a * fb * fc / ((fa - fb) * (fa - fc))
        + b * fa * fc / ((fb - fa) * (fb - fc))
        + c * fa * fb / ((fc - fa) * (fc - fb));
    } else {
      s = b - fb * (b - a) / (fb - fa);
    }

    const lo = (3 * a + b) / 4;
    const outOfRange = !(((s > lo) && (s < b)) || ((s < lo) && (s > b)));
    const cond2 = mflag && Math.abs(s - b) >= Math.abs(b - c) / 2;
    const cond3 = !mflag && Math.abs(s - b) >= Math.abs(c - d) / 2;
    const cond4 = mflag && Math.abs(b - c) < tol;
    const cond5 = !mflag && Math.abs(c - d) < tol;

    if (outOfRange || cond2 || cond3 || cond4 || cond5) {
      s = (a + b) / 2;
      mflag = true;
    } else {
      mflag = false;
    }

    const fs = Z(s);
    d = c; c = b; fc = fb;
    if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
    if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
  }
  return b;
}

function computeCurveAndZeros(tStart, tEnd, samples = 600) {
  const step = (tEnd - tStart) / samples;
  const points = [];
  const zeros = [];

  let prevT = tStart;
  let prevZ = Z(tStart);
  points.push({ t: prevT, z: prevZ });

  for (let i = 1; i <= samples; i++) {
    const t = tStart + i * step;
    const z = Z(t);
    points.push({ t, z });
    if (prevZ * z < 0) {
      zeros.push(brentZero(prevT, t));
    }
    prevT = t;
    prevZ = z;
  }
  return { points, zeros };
}

// =========================================================================
// Canvas rendering
// =========================================================================

const canvas = document.getElementById('zetaCanvas');
const ctx = canvas.getContext('2d');

function resizeCanvas() {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = rect.width * dpr;
  canvas.height = rect.height * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function drawPlot(points, zeros, tStart, tEnd) {
  const rect = canvas.getBoundingClientRect();
  const w = rect.width;
  const h = rect.height;
  const padL = 44, padR = 12, padT = 16, padB = 26;
  const plotW = w - padL - padR;
  const plotH = h - padT - padB;

  ctx.clearRect(0, 0, w, h);

  let zMin = Infinity, zMax = -Infinity;
  for (const p of points) {
    if (p.z < zMin) zMin = p.z;
    if (p.z > zMax) zMax = p.z;
  }
  const zPad = (zMax - zMin) * 0.1 || 1;
  zMin -= zPad; zMax += zPad;

  const xOf = (t) => padL + ((t - tStart) / (tEnd - tStart)) * plotW;
  const yOf = (z) => padT + (1 - (z - zMin) / (zMax - zMin)) * plotH;

  // grid
  ctx.strokeStyle = '#1B2530';
  ctx.lineWidth = 1;
  ctx.beginPath();
  const zeroY = yOf(0);
  ctx.moveTo(padL, zeroY);
  ctx.lineTo(padL + plotW, zeroY);
  ctx.stroke();

  // axis labels
  ctx.fillStyle = '#8B98A6';
  ctx.font = '11px "IBM Plex Mono", monospace';
  ctx.fillText(tStart.toFixed(0), padL, h - 8);
  ctx.fillText(tEnd.toFixed(0), padL + plotW - 24, h - 8);
  ctx.fillText('Z(t)=0', 4, zeroY + 4);

  // curve
  ctx.strokeStyle = '#E8B34C';
  ctx.lineWidth = 1.6;
  ctx.beginPath();
  points.forEach((p, i) => {
    const x = xOf(p.t), y = yOf(p.z);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  });
  ctx.stroke();

  // zero markers
  ctx.fillStyle = '#4CC9C0';
  for (const zt of zeros) {
    const x = xOf(zt), y = yOf(0);
    ctx.beginPath();
    ctx.arc(x, y, 3.5, 0, Math.PI * 2);
    ctx.fill();
  }
}

// =========================================================================
// UI wiring
// =========================================================================

const tStartInput = document.getElementById('tStart');
const tEndInput = document.getElementById('tEnd');
const tStartVal = document.getElementById('tStartVal');
const tEndVal = document.getElementById('tEndVal');
const zeroCountEl = document.getElementById('zeroCount');
const computeTimeEl = document.getElementById('computeTime');
const terminalEl = document.getElementById('terminal');

function render() {
  let tStart = parseFloat(tStartInput.value);
  let tEnd = parseFloat(tEndInput.value);
  if (tEnd <= tStart + 5) {
    tEnd = tStart + 5;
    tEndInput.value = tEnd;
  }

  tStartVal.textContent = tStart.toFixed(0);
  tEndVal.textContent = tEnd.toFixed(0);

  const t0 = performance.now();
  const { points, zeros } = computeCurveAndZeros(tStart, tEnd);
  const elapsed = performance.now() - t0;

  drawPlot(points, zeros, tStart, tEnd);

  zeroCountEl.textContent = zeros.length;
  computeTimeEl.textContent = `${elapsed.toFixed(2)} ms (brauzerde)`;

  const lines = [
    '=========================================================',
    ' RIEMANN ZETA ENGINE — JAVASCRIPT JANLY DEMO',
    '=========================================================',
    `[Range] t in [${tStart.toFixed(1)}, ${tEnd.toFixed(1)}]`,
    `[Precision] double (64-bit) — brauzer demosy üçin`,
    '',
  ];
  zeros.forEach((z, i) => {
    lines.push(`Zero #${i + 1} -> s = 0.5 + ${z.toFixed(10)}i`);
  });
  lines.push('');
  lines.push('---------------------------------------------------------');
  lines.push(`Total zeros found: ${zeros.length}`);
  lines.push(`Elapsed time: ${elapsed.toFixed(3)} ms`);
  lines.push('=========================================================');

  terminalEl.textContent = lines.join('\n');
}

let pending = false;
function scheduleRender() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => {
    pending = false;
    render();
  });
}

tStartInput.addEventListener('input', scheduleRender);
tEndInput.addEventListener('input', scheduleRender);
window.addEventListener('resize', () => { resizeCanvas(); scheduleRender(); });

resizeCanvas();
render();
