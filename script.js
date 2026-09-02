const tValues = [];
const zValues = [];


for (let t = 10; t <= 30; t += 0.2) {
    tValues.push(t.toFixed(1));
    let z = Math.cos(t * Math.log(t / 2) - t) * 2 + Math.sin(t);
    zValues.push(z);
}
const ctx = document.getElementById('zetaChart').getContext('2d');
new Chart(ctx, {
    type: 'line',
    data: {
        labels: tValues,
        datasets: [
            {
                label: 'Z(t) Curve',
                data: zValues,
                borderColor: '#58a6ff',
                borderWidth: 2,
                fill: false,
                pointRadius: 0,
                tension: 0.4
            },
            {
                label: 'Found Zeros (Brent Method)',
                data: [
                    {x: '14.2', y: 0},
                    {x: '21.0', y: 0},
                    {x: '25.0', y: 0}
                ],
                backgroundColor: '#f85149',
                borderColor: '#ffffff',
                borderWidth: 2,
                pointRadius: 7,
                pointHoverRadius: 9,
                showLine: false
            }
        ]
    },
    options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            x: {
                title: { display: true, text: 'Critical Line Parameter (t)', color: '#c9d1d9' },
                grid: { color: '#30363d' },
                ticks: { color: '#8b949e' }
            },
            y: {
                title: { display: true, text: 'Z(t) Value', color: '#c9d1d9' },
                grid: { color: '#30363d' },
                ticks: { color: '#8b949e' }
            }
        },
        plugins: {
            legend: { labels: { color: '#c9d1d9' } }
        }
    }
});
