/* =========================================================
   Numerical routines
   =========================================================
   - Regularized lower incomplete gamma P(a,x)
   - Central chi-square quantile
   - Noncentral chi-square CDF via Poisson mixture
   - Root solve for noncentrality λ giving requested power

   This keeps the calculator self-contained:
   no JavaScript library is required.
   ========================================================= */


/* ---------------------------------------------------------
   Log Gamma Function
   --------------------------------------------------------- */

function logGamma(z) {

    const c = [
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ];

    if (z < 0.5) {
        return (
            Math.log(Math.PI) -
            Math.log(Math.sin(Math.PI * z)) -
            logGamma(1 - z)
        );
    }

    z -= 1;

    let x = 0.99999999999980993;

    for (let i = 0; i < c.length; i++) {
        x += c[i] / (z + i + 1);
    }

    const t = z + c.length - 0.5;

    return (
        0.5 * Math.log(2 * Math.PI) +
        (z + 0.5) * Math.log(t) -
        t +
        Math.log(x)
    );
}


/* ---------------------------------------------------------
   Regularized Lower Incomplete Gamma Function
   P(a,x)
   --------------------------------------------------------- */

function gammaP(a, x) {

    if (x <= 0) {
        return 0;
    }

    /* Series representation */
    if (x < a + 1) {

        let ap = a;
        let sum = 1 / a;
        let del = sum;

        for (let n = 1; n <= 10000; n++) {

            ap += 1;
            del *= x / ap;
            sum += del;

            if (
                Math.abs(del) <
                Math.abs(sum) * 3e-14
            ) {
                break;
            }
        }

        return (
            sum *
            Math.exp(
                -x +
                a * Math.log(x) -
                logGamma(a)
            )
        );
    }


    /* Continued fraction representation */

    let b = x + 1 - a;
    let c = 1 / 1e-300;
    let d = 1 / b;
    let h = d;

    for (let i = 1; i <= 10000; i++) {

        const an = -i * (i - a);
        const bb = b + 2 * i;

        d = an * d + bb;

        if (Math.abs(d) < 1e-300) {
            d = 1e-300;
        }

        c = bb + an / c;

        if (Math.abs(c) < 1e-300) {
            c = 1e-300;
        }

        d = 1 / d;

        const del = d * c;

        h *= del;

        if (Math.abs(del - 1) < 3e-14) {
            break;
        }
    }

    return (
        1 -
        Math.exp(
            -x +
            a * Math.log(x) -
            logGamma(a)
        ) * h
    );
}


/* ---------------------------------------------------------
   Chi-Square CDF
   --------------------------------------------------------- */

function chiSquareCDF(x, df) {

    return gammaP(
        df / 2,
        x / 2
    );
}


/* ---------------------------------------------------------
   Normal Inverse CDF
   Acklam Approximation
   --------------------------------------------------------- */

function normalInv(p) {

    const a = [
        -39.6968302866538,
        220.946098424521,
        -275.928510446969,
        138.357751867269,
        -30.6647980661472,
        2.50662827745924
    ];

    const b = [
        -54.4760987982241,
        161.585836858041,
        -155.698979859887,
        66.8013118877197,
        -13.2806815528857
    ];

    const c = [
        -0.00778489400243029,
        -0.322396458041136,
        -2.40075827716184,
        -2.54973253934373,
        4.37466414146497,
        2.93816398269878
    ];

    const d = [
        0.00778469570904146,
        0.32246712907004,
        2.445134137143,
        3.75440866190742
    ];

    const pl = 0.02425;
    const ph = 1 - pl;

    if (p <= 0) {
        return -Infinity;
    }

    if (p >= 1) {
        return Infinity;
    }

    let q;
    let r;


    /* Lower region */

    if (p < pl) {

        q = Math.sqrt(-2 * Math.log(p));

        return (
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        );
    }


    /* Upper region */

    if (p > ph) {

        q = Math.sqrt(-2 * Math.log(1 - p));

        return -(
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        );
    }


    /* Central region */

    q = p - 0.5;
    r = q * q;

    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q) /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    );
}


/* ---------------------------------------------------------
   Chi-Square Quantile
   --------------------------------------------------------- */

function chiSquareQuantile(p, df) {

    let lo = 0;
    let hi = Math.max(df, 1);

    while (chiSquareCDF(hi, df) < p) {
        hi *= 2;
    }

    for (let i = 0; i < 100; i++) {

        const mid = (lo + hi) / 2;

        if (chiSquareCDF(mid, df) < p) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return (lo + hi) / 2;
}


/* ---------------------------------------------------------
   Noncentral Chi-Square CDF
   --------------------------------------------------------- */

function noncentralChiCDF(x, df, ncp) {

    if (ncp < 1e-10) {
        return chiSquareCDF(x, df);
    }

    const mu = ncp / 2;
    const k0 = Math.floor(mu);

    let weight =
        Math.exp(
            -mu +
            k0 * Math.log(mu) -
            logGamma(k0 + 1)
        );

    let sum = 0;

    /* Start at the mode */

    sum +=
        weight *
        gammaP(
            df / 2 + k0,
            x / 2
        );


    /* Move upward */

    let w = weight;

    for (
        let k = k0 + 1;
        k < k0 + 2000;
        k++
    ) {

        w *= mu / k;

        const term =
            w *
            gammaP(
                df / 2 + k,
                x / 2
            );

        sum += term;

        if (
            k > mu + 30 * Math.sqrt(mu + 1) &&
            term < 1e-14
        ) {
            break;
        }
    }


    /* Move downward */

    w = weight;

    for (
        let k = k0 - 1;
        k >= 0;
        k--
    ) {

        w *= (k + 1) / mu;

        const term =
            w *
            gammaP(
                df / 2 + k,
                x / 2
            );

        sum += term;

        if (
            k < mu - 30 * Math.sqrt(mu + 1) &&
            term < 1e-14
        ) {
            break;
        }
    }

    return Math.max(
        0,
        Math.min(1, sum)
    );
}


/* ---------------------------------------------------------
   Required Noncentrality Parameter
   --------------------------------------------------------- */

function requiredNcp(alpha, power, df) {

    /* Critical value under the null */
    const critical =
        chiSquareQuantile(
            1 - alpha,
            df
        );

    /* Desired CDF */
    const targetCDF = 1 - power;

    let lo = 0;
    let hi = Math.max(10, df);


    /* Find an upper bound */

    while (
        noncentralChiCDF(
            critical,
            df,
            hi
        ) > targetCDF
    ) {
        hi *= 2;
    }


    /* Binary search */

    for (let i = 0; i < 80; i++) {

        const mid = (lo + hi) / 2;

        if (
            noncentralChiCDF(
                critical,
                df,
                mid
            ) > targetCDF
        ) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return (lo + hi) / 2;
}


/* ---------------------------------------------------------
   Model Degrees of Freedom
   --------------------------------------------------------- */

function modelDf(p, k) {

    const b =
        p * (p + 1) / 2;

    const ncor =
        k * (k - 1) / 2;

    const a =
        (p - k) +
        p +
        k +
        ncor;

    return Math.round(b - a);
}


/* ---------------------------------------------------------
   RMSEA Sample Size Calculation
   --------------------------------------------------------- */

function calculate() {

    const rmsea =
        Number(
            document.getElementById('rmsea').value
        );

    const alpha =
        Number(
            document.getElementById('alpha').value
        );

    const power =
        Number(
            document.getElementById('power').value
        ) / 100;

    const df =
        Math.round(
            Number(
                document.getElementById('df').value
            )
        );

    const drop =
        Number(
            document.getElementById('dropout').value
        ) / 100;

    const warning =
        document.getElementById('warning');

    warning.style.display = 'none';


    /* Validate inputs */

    if (
        !(rmsea > 0) ||
        !(alpha > 0 && alpha < 1) ||
        !(power > 0 && power < 1) ||
        !(df > 0) ||
        !(drop >= 0 && drop < 1)
    ) {

        warning.textContent =
            'Please enter valid values. RMSEA and df must be positive; ' +
            'α and power must be between 0 and 1; dropout must be below 100%.';

        warning.style.display = 'block';

        return;
    }


    try {

        const ncp =
            requiredNcp(
                alpha,
                power,
                df
            );

        const n =
            ncp /
            (rmsea * rmsea * df) +
            1;

        const nd =
            n /
            (1 - drop);


        document.getElementById('rDf')
            .textContent = df;

        document.getElementById('rNcp')
            .textContent = ncp.toFixed(3);

        document.getElementById('rN')
            .textContent = Math.ceil(n);

        document.getElementById('rDrop')
            .textContent = Math.ceil(nd);

    } catch (e) {

        warning.textContent =
            'The numerical calculation could not converge for these inputs. ' +
            'Try conventional values such as RMSEA 0.05, α 0.05, and power 80%.';

        warning.style.display = 'block';
    }
}


/* ---------------------------------------------------------
   SEM Degrees of Freedom Calculation
   --------------------------------------------------------- */

function calculateDf() {

    const p =
        Math.round(
            Number(
                document.getElementById('items').value
            )
        );

    const vals = [
        'loadings',
        'errors',
        'exoVars',
        'exoCov',
        'paths',
        'disturbances'
    ].map(
        id =>
            Math.round(
                Number(
                    document.getElementById(id).value
                )
            )
    );

    const warning =
        document.getElementById('warning');


    /* Validate input */

    if (
        !(p >= 2) ||
        vals.some(v => !(v >= 0))
    ) {

        warning.textContent =
            'Please enter valid non-negative parameter counts ' +
            'and at least 2 observed indicators.';

        warning.style.display = 'block';

        return;
    }


    /* Number of unique covariance elements */

    const unique =
        p * (p + 1) / 2;


    /* Number of free parameters */

    const free =
        vals.reduce(
            (a, b) => a + b,
            0
        );


    /* Degrees of freedom */

    const df =
        Math.round(
            unique - free
        );


    document.getElementById('uniqueResult')
        .textContent = unique;

    document.getElementById('freeResult')
        .textContent = free;

    document.getElementById('dfResult')
        .textContent = df;

    document.getElementById('df')
        .value = df;


    /* Check whether model is identified */

    if (df <= 0) {

        warning.textContent =
            'This specification produces non-positive degrees of freedom.';

        warning.style.display = 'block';

    } else {

        warning.style.display = 'none';
    }
}


/* ---------------------------------------------------------
   Event Listeners
   --------------------------------------------------------- */

document
    .getElementById('calculate')
    .addEventListener(
        'click',
        calculate
    );


document
    .getElementById('calcDf')
    .addEventListener(
        'click',
        calculateDf
    );


document
    .getElementById('reset')
    .addEventListener(
        'click',
        () => {

            setTimeout(() => {

                document.getElementById('dfResult')
                    .textContent = '241';

                document.getElementById('rDf')
                    .textContent = '241';

                document.getElementById('rNcp')
                    .textContent = '—';

                document.getElementById('rN')
                    .textContent = '—';

                document.getElementById('rDrop')
                    .textContent = '—';

                document.getElementById('warning')
                    .style.display = 'none';

            }, 0);
        }
    );


/* ---------------------------------------------------------
   Initial Calculation
   --------------------------------------------------------- */

calculate();
