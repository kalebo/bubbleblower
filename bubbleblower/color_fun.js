// Derived from Martin Ankerl's blog about random programmatic colors.
function map_color(i){
    phi_conj = 0.618033988749895;
    h = ((i * .944) + phi_conj) % 1;
    s = ((i * .7) + phi_conj) % 1;
    v = ((i * .1) + phi_conj) % 1;
    sa = (s > .5) ? s : s + .4;
    va = (v > .5) ? v : v + .4;
    return hvs2rgb(h, va ,sa);
}

function hvs2rgb(h, v, s){
    var h_sextuple = h * 6;
    var h_int = Math.floor(h_sextuple);

    var f = h_sextuple - h_int;
    var p = v * (1 - s);
    var q = v * (1 - f * s);
    var t = v * (1 - (1 - f) * s);

    var r, g, b;
    switch (h_int) {
    case 0: [r, g, b] = [v, t, p];
        break;
    case 1: [r, g, b] = [q, v, p];
        break;
    case 2: [r, g, b] = [p, v, t];
        break;
    case 3: [r, g, b] = [p, q, v];
        break;
    case 4: [r, g, b] = [t, p, v];
        break;
    case 5: [r, g, b] = [v, p, q];
        break;
    }

    rgb = [r,g,b]
        .map(function(x){return x * 256;})
        .map(Math.floor)
        .map(function(x){return x.toString(16);})
        .map(function(x){return (x.length < 2) ? '0' + x : x;})
        .join("");

    return rgb;



}

function change_bgcolor(rgb) {
    document.body.style.backgroundColor = rgb;
}


