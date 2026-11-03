"""The single HTML page the web front end serves.

Kept in its own module so `webapp.py` stays about HTTP and the model, and
this stays about markup.  It is one string on purpose: no build step, no
bundler, no CDN, nothing to install - which is the same rule the rest of
the project follows.

The JavaScript here does three things only: collect the form into JSON,
POST it, and drop the returned HTML and SVG into the page.  Every number
on screen was computed by the Python model.
"""
from __future__ import annotations

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DBF 2026-27 optimiser</title>
<style>
  :root{
    --bg:#f4f6f9; --panel:#fff; --line:#dde3ea; --ink:#1c2430;
    --muted:#6b7580; --accent:#2f6fdb; --good:#2e9e5b; --bad:#d9534f;
  }
  *{box-sizing:border-box}
  body{margin:0;font:14px/1.45 system-ui,"Segoe UI",sans-serif;
       background:var(--bg);color:var(--ink)}
  header{background:var(--panel);border-bottom:1px solid var(--line);
         padding:10px 16px;display:flex;align-items:center;gap:16px;
         position:sticky;top:0;z-index:10;flex-wrap:wrap}
  header h1{font-size:15px;margin:0;font-weight:600}
  #status{margin-left:auto;font-size:13px;color:var(--muted)}
  #status .ok{color:var(--good);font-weight:600}
  #status .bad{color:var(--bad);font-weight:600}
  main{display:flex;gap:14px;padding:14px;align-items:flex-start}
  #left{flex:0 0 380px;position:sticky;top:60px}
  #right{flex:1;min-width:0}
  .card{background:var(--panel);border:1px solid var(--line);
        border-radius:8px;margin-bottom:12px;overflow:hidden}
  .card h2{font-size:12px;text-transform:uppercase;letter-spacing:.06em;
           margin:0;padding:9px 12px;border-bottom:1px solid var(--line);
           color:var(--muted);background:#fafbfc}
  .card .body{padding:10px 12px}
  .row{display:flex;align-items:center;gap:8px;margin:4px 0}
  .row label{flex:0 0 150px;font-size:13px}
  .row input,.row select{flex:1;min-width:0;padding:4px 6px;font:13px inherit;
    border:1px solid var(--line);border-radius:4px;background:#fff}
  .row .unit{flex:0 0 74px;font-size:11px;color:var(--muted)}
  .note{font-size:12px;color:var(--accent);margin:6px 0 2px}
  .hint{font-size:12px;color:var(--muted);margin:6px 0 0}
  button{padding:6px 12px;font:13px inherit;border:1px solid var(--line);
         border-radius:5px;background:#fff;cursor:pointer}
  button.primary{background:var(--accent);color:#fff;border-color:var(--accent)}
  button:disabled{opacity:.5;cursor:default}
  .tabs{display:flex;gap:2px;border-bottom:1px solid var(--line);
        background:#fafbfc;flex-wrap:wrap}
  .tabs button{border:0;border-radius:0;background:transparent;
               padding:9px 13px;color:var(--muted);font-size:13px}
  .tabs button.on{color:var(--ink);font-weight:600;
                  box-shadow:inset 0 -2px 0 var(--accent);background:#fff}
  .pane{display:none;padding:12px}
  .pane.on{display:block}
  pre{margin:0;font:12px/1.4 Consolas,ui-monospace,monospace;
      white-space:pre;overflow:auto;max-height:74vh}
  table{border-collapse:collapse;width:100%;font-size:13px}
  th,td{text-align:left;padding:4px 8px;border-bottom:1px solid #eef1f5}
  th{color:var(--muted);font-weight:600;font-size:12px}
  td.num{font-family:Consolas,monospace}
  .block{border:1px solid var(--line);border-radius:6px;margin-bottom:12px}
  .block h3{margin:0;padding:8px 12px;font-size:13px;background:#fafbfc;
            border-bottom:1px solid var(--line)}
  .block .formula{padding:8px 12px;font:13px Consolas,monospace;
                  background:#f7f9fc;color:#24425e}
  .block .result{padding:8px 12px;font:600 14px Consolas,monospace;
                 border-top:1px solid var(--line)}
  .block .why{padding:0 12px 10px;font-size:12px;color:var(--muted)}
  .sat{color:var(--bad);font-weight:600}
  .err{color:var(--bad)}
  details{border-top:1px solid var(--line)}
  details summary{padding:7px 12px;cursor:pointer;font-size:12px;
                  color:var(--muted);background:#fafbfc}
  .adv{padding:8px 12px}
  svg{max-width:100%;height:auto;display:block}
  @media (max-width:900px){
    main{flex-direction:column}
    #left{flex:1 1 auto;position:static;width:100%}
  }
</style>
</head>
<body>
<header>
  <h1>DBF 2026-27 optimiser</h1>
  <button class="primary" id="go">Evaluate</button>
  <select id="mode"><option>quick</option><option>default</option>
    <option>thorough</option></select>
  <button id="search">Run search</button>
  <label style="font-size:13px"><input type="checkbox" id="imperial" checked>
    imperial</label>
  <div id="status">loading…</div>
</header>

<main>
<div id="left">
  <div class="card"><h2>Wing</h2><div class="body" id="panel-wing"></div></div>
  <div class="card"><h2>Propulsion</h2><div class="body" id="panel-prop"></div>
    <div class="body" style="padding-top:0"><div class="note" id="packnote">
    </div></div></div>
  <div class="card"><h2>Payload</h2><div class="body" id="panel-payload"></div></div>
  <div class="card"><h2>Conditions on the day</h2>
    <div class="body" id="panel-env"></div></div>
  <div class="card"><h2>Search space</h2>
    <div class="body hint" style="color:#a0522d">The ranges <b>Run search</b>
      explores. The design fields above are one aeroplane that
      <b>Evaluate</b> scores exactly &mdash; a search ignores them and
      writes its winner back into them.</div>
    <div id="space"></div></div>
  <div class="card"><h2>Advanced</h2>
    <div class="body hint">Every remaining number in Rules, Course and the
      build standard. Generated from the dataclasses, so nothing is
      unreachable.</div>
    <div id="advanced"></div></div>
</div>

<div id="right">
  <div class="card">
    <div class="tabs" id="tabs"></div>
    <div class="pane on" id="p-score"></div>
    <div class="pane" id="p-report"><pre id="report"></pre></div>
    <div class="pane" id="p-dims"><table id="dims"></table></div>
    <div class="pane" id="p-charts">
      <div class="row"><label>chart</label><select id="chart"></select>
        <span class="unit"></span></div>
      <div id="chartbox"></div></div>
    <div class="pane" id="p-3v"><div id="tvbox"></div></div>
    <div class="pane" id="p-bom"><pre id="bom"></pre></div>
    <div class="pane" id="p-search"><pre id="searchlog">Not run yet.</pre></div>
    <div class="pane" id="p-rules"><pre id="rules">loading…</pre></div>
    <div class="pane" id="p-cat">
      <div class="row"><label>filter</label>
        <input id="catq" placeholder="e.g. trim, spar, pack"></div>
      <pre id="cat"></pre></div>
  </div>
</div>
</main>

<script>
const $ = s => document.querySelector(s);
let META = null, BUSY = false;

const TABS = [["score","Score"],["report","Report"],["dims","Dimensions"],
              ["charts","Charts"],["3v","Three-view"],["bom","Parts list"],
              ["search","Search"],["rules","Rules"],
              ["cat","Function catalogue"]];

function buildTabs(){
  $("#tabs").innerHTML = TABS.map(([id,label],i)=>
    `<button data-t="${id}" class="${i===0?'on':''}">${label}</button>`).join("");
  $("#tabs").onclick = e => {
    const b = e.target.closest("button"); if(!b) return;
    document.querySelectorAll("#tabs button").forEach(x=>x.classList.remove("on"));
    document.querySelectorAll(".pane").forEach(x=>x.classList.remove("on"));
    b.classList.add("on");
    $("#p-"+b.dataset.t).classList.add("on");
  };
}

function row(id, label, value, unit, opts){
  if(opts){
    const o = opts.map(v => {
      const [k, t] = Array.isArray(v) ? v : [v, v];
      return `<option value="${k}"${k===value?" selected":""}>${t}</option>`;
    }).join("");
    return `<div class="row"><label for="${id}">${label}</label>
      <select id="${id}">${o}</select><span class="unit">${unit||""}</span></div>`;
  }
  return `<div class="row"><label for="${id}">${label}</label>
    <input id="${id}" value="${value}"><span class="unit">${unit||""}</span></div>`;
}

function buildForm(m){
  // The form opens on the same design as the desktop GUI - the server's
  // opening design, read from examples/optimised_v1.json - rather than
  // values typed in here, which went stale as the optimum moved.
  const o = m.opening;
  $("#panel-wing").innerHTML =
    '<p class="hint">opens on: ' + esc(o.source) + '</p>' +
    row("airfoil","airfoil",o.airfoil,"",m.airfoils) +
    row("span_m","span",(+o.span_m).toFixed(4),"m") +
    row("aspect_ratio","aspect ratio",(+o.aspect_ratio).toFixed(4),
        "= span/chord") +
    row("taper","taper ratio",(+o.taper).toFixed(4),"tip/root") +
    row("winglet_h_m","winglet height",(+o.winglet_h_m).toFixed(4),
        "m each tip, 0 none") +
    row("tail_vol_h","tail volume Vh",(+o.tail_vol_h).toFixed(3),
        "0 = just its requirements");
  $("#panel-prop").innerHTML =
    row("motor","motor",o.motor,"",m.motors) +
    row("prop","propeller",o.prop,"",m.props) +
    row("n_motors","motor count",o.n_motors,"max 2") +
    row("pack","battery pack",o.pack,"",
        m.packs.map(p=>[p.key,p.label])) +
    row("n_systems","propulsion systems",o.n_systems,"one pack each") +
    row("wing_prop_layout","wing props",o.wing_prop_layout,
        "two or more motors",["mid","tip","stagger"]);
  $("#panel-payload").innerHTML =
    row("sensor_mass_kg","sensor mass (declared max)",
        (+o.sensor_mass_kg).toFixed(4),"kg") +
    row("n_simulators","simulators",o.n_simulators,"each = loaded box") +
    row("sensor_m3_kg","M3 sensor mass",
        (+(o.sensor_m3_kg || o.sensor_mass_kg)).toFixed(4),"kg, <= max") +
    row("drop_height_in","drop height",o.drop_height_in,"in, max 60");
  $("#panel-env").innerHTML =
    '<p class="hint">' + esc(m.site) + '</p>' +
    row("conditions","preset",m.default_conditions,"",
        Object.keys(m.conditions)) +
    row("wind_mps","wind",m.conditions[m.default_conditions][2],"m/s") +
    row("altitude_m","field elevation",
        m.conditions[m.default_conditions][0],"m") +
    row("temp_c","temperature",
        m.conditions[m.default_conditions][1],"deg C") +
    row("rho","air density","1.225","kg/m3") +
    '<p class="note" id="rhonote"></p>';
  const sp = m.space;
  $("#space").innerHTML = '<div class="adv">' +
    row("sp:span_lo","span min",sp.span[0],"m") +
    row("sp:span_hi","span max",sp.span[1],"capped by the 6 ft rule") +
    row("sp:ar_lo","aspect ratio min",sp.ar[0],"") +
    row("sp:ar_hi","aspect ratio max",sp.ar[1],"") +
    row("sp:taper_lo","taper min",sp.taper[0],"") +
    row("sp:taper_hi","taper max",sp.taper[1],"") +
    row("sp:sensor_lo","sensor min",sp.sensor[0],"kg") +
    row("sp:sensor_hi","sensor max",sp.sensor[1],"kg") +
    row("sp:m3_frac_lo","M3 fraction min",sp.m3_frac[0],"of max") +
    row("sp:m3_frac_hi","M3 fraction max",sp.m3_frac[1],"") +
    row("sp:drop_lo","drop height min",sp.drop[0],"in") +
    row("sp:drop_hi","drop height max",sp.drop[1],"max 60") +
    row("sp:motor_counts","motor counts",sp.motor_counts,"max 4 (team cap)") +
    row("sp:system_counts","propulsion systems",sp.system_counts,"") +
    '<p class="hint">Airfoils (' + sp.catalogue[0] + '), motors (' +
    sp.catalogue[1] + '), propellers (' + sp.catalogue[2] +
    ') and packs (' + sp.catalogue[3] +
    ') are searched in full from the catalogues.</p></div>';
  $("#advanced").innerHTML = m.advanced.map(g=>
    `<details><summary>${g.title} (${g.fields.length})</summary>
     <div class="adv">${g.fields.map(f=>
       row("adv:"+g.owner+":"+f.name, f.name, f.value, "")).join("")}</div>
     </details>`).join("");
  $("#chart").innerHTML = m.charts.map(c=>`<option>${c}</option>`).join("");
  $("#pack").onchange = packNote;
  $("#chart").onchange = evaluate;
  $("#conditions").onchange = () => {
    const c = META.conditions[$("#conditions").value];
    if(!c) return;
    $("#altitude_m").value = c[0];
    $("#temp_c").value = c[1];
    $("#wind_mps").value = c[2];
    $("#rho").value = "";        // blank: let the server re-derive it
    evaluate();
  };
  packNote();
}

function packNote(){
  const k = $("#pack").value;
  const p = META.packs.find(x=>x.key===k);
  $("#packnote").textContent = p ? p.label : "";
}

function collect(){
  const d = {};
  // If an id here ever collides with a container's id, $("#x") returns the
  // container, .value is undefined, the key is dropped from the JSON and
  // the server silently substitutes its own default - which once meant a
  // different propeller than the one on screen.  Fail loudly instead.
  ["airfoil","span_m","aspect_ratio","taper","motor","prop","n_motors",
   "pack","n_systems","sensor_mass_kg","n_simulators","sensor_m3_kg",
   "drop_height_in","winglet_h_m","tail_vol_h","wing_prop_layout"]
    .forEach(k => {
      const el = $("#"+k);
      if(!el || el.value === undefined)
        throw new Error("form field '"+k+"' is missing or shadowed");
      d[k] = el.value;
    });
  const env = {};
  ["wind_mps","altitude_m","temp_c","rho"].forEach(k => env[k] = $("#"+k).value);
  const out = {design:d, env:env, rules:{}, course:{}, build:{}, space:{},
               imperial:$("#imperial").checked, chart:$("#chart").value};
  document.querySelectorAll('[id^="sp:"]').forEach(el=>{
    out.space[el.id.slice(3)] = el.value;
  });
  document.querySelectorAll('[id^="adv:"]').forEach(el=>{
    const [,owner,name] = el.id.split(":");
    out[owner][name] = el.value;
  });
  return out;
}

function esc(s){ return String(s).replace(/[&<>]/g,
  c => ({"&":"&amp;","<":"&lt;",">":"&gt;"}[c])); }

function renderScore(blocks){
  $("#p-score").innerHTML = blocks.map(b=>`
    <div class="block">
      <h3>${esc(b.title)}</h3>
      <div class="formula">${esc(b.formula)}</div>
      <table>${b.rows.map(r=>
        `<tr><td>${esc(r[0])}</td><td class="num">${esc(r[1])}</td>
         <td style="color:var(--muted);font-size:12px">${esc(r[2])}</td></tr>`
        ).join("")}</table>
      <div class="result">${esc(b.result)}</div>
      ${b.note ? `<div class="why${b.note.startsWith("SATURATED")?" sat":""}">
         ${esc(b.note)}</div>` : ""}
    </div>`).join("");
}

function renderDims(rows){
  $("#dims").innerHTML =
    "<tr><th>Quantity</th><th>SI</th><th>Imperial</th></tr>" +
    rows.map(r => r[0] === "" ? '<tr><td colspan="3">&nbsp;</td></tr>'
      : `<tr><td>${esc(r[0])}</td><td class="num">${esc(r[1])}</td>
         <td class="num">${esc(r[2])}</td></tr>`).join("");
}

async function evaluate(){
  if(BUSY) return;
  BUSY = true; $("#go").disabled = true;
  $("#status").textContent = "evaluating…";
  try{
    const r = await fetch("/api/evaluate", {method:"POST",
      headers:{"Content-Type":"application/json"},
      body: JSON.stringify(collect())});
    const j = await r.json();
    // Clear anything the response does not carry. Guarding each panel
    // with `if (j.x)` leaves the PREVIOUS design's dimensions, charts,
    // three-view and score on screen when a request fails - a complete,
    // plausible aeroplane that is not the one in the form, with nothing
    // to say so.
    const gone = j.ok === false && !j.dimensions
      ? "This design could not be evaluated." : "";
    $("#report").textContent = j.report || gone;
    $("#bom").textContent = j.bom || gone;
    if(j.score_blocks) renderScore(j.score_blocks);
    else $("#p-score").innerHTML = '<p class="err">' + esc(gone) + '</p>';
    if(j.dimensions) renderDims(j.dimensions);
    else $("#dims").innerHTML = '<tr><td class="err">' + esc(gone) +
      '</td></tr>';
    $("#chartbox").innerHTML = j.chart_svg ||
      ('<p class="err">' + esc(gone) + '</p>');
    $("#tvbox").innerHTML = j.three_view_svg ||
      ('<p class="err">' + esc(gone) + '</p>');
    if(j.env) $("#rhonote").textContent =
      "= " + j.env.density_altitude_ft.toFixed(0) + " ft density altitude, "
      + j.env.isa_pct.toFixed(1) + "% of sea-level ISA";
    const s = j.summary;
    if(s){
      $("#status").innerHTML =
        `<span class="${j.ok?"ok":"bad"}">${j.ok?"OK":"INFEASIBLE"}</span>` +
        ` &nbsp; M2 ${s.m2.toFixed(3)} &nbsp; M3 ${s.m3.toFixed(3)}` +
        ` &nbsp; total ${s.total.toFixed(3)} &nbsp;|&nbsp; ` +
        `${s.mtow_lb.toFixed(1)} lb, ${s.watt_hours.toFixed(1)} Wh, ` +
        (j.structure_ok ? "structure OK" : "SPAR WILL NOT CLOSE");
    } else {
      $("#status").innerHTML = '<span class="bad">failed</span>';
    }
    if(!j.ok && j.failures && j.failures.length){
      $("#report").textContent =
        "CONSTRAINTS VIOLATED:\n  - " + j.failures.join("\n  - ") +
        "\n\n" + ($("#report").textContent || "");
    }
  }catch(e){
    $("#status").innerHTML = '<span class="bad">' + esc(e) + '</span>';
  }finally{
    BUSY = false; $("#go").disabled = false;
  }
}

async function startSearch(){
  const body = collect();
  body.mode = $("#mode").value;
  $("#search").disabled = true;
  document.querySelector('[data-t="search"]').click();
  $("#searchlog").textContent = "starting…";
  await fetch("/api/search", {method:"POST",
    headers:{"Content-Type":"application/json"}, body: JSON.stringify(body)});
  poll();
}

async function poll(){
  const j = await (await fetch("/api/search")).json();
  $("#searchlog").textContent = j.log || "";
  if(j.state === "running"){ setTimeout(poll, 1500); return; }
  $("#search").disabled = false;
  if(j.state === "done" && j.best){
    const b = j.best;
    $("#airfoil").value = b.airfoil; $("#motor").value = b.motor;
    $("#prop").value = b.prop; $("#pack").value = b.pack;
    $("#n_motors").value = b.n_motors; $("#n_systems").value = b.n_systems;
    $("#span_m").value = (+b.span_m).toFixed(4);
    $("#aspect_ratio").value = (b.span_m / b.chord_m).toFixed(4);
    $("#taper").value = (+b.taper).toFixed(4);
    $("#sensor_mass_kg").value = (+b.sensor_mass_kg).toFixed(4);
    $("#n_simulators").value = b.n_simulators;
    $("#sensor_m3_kg").value =
      (+(b.sensor_m3_kg || b.sensor_mass_kg)).toFixed(4);
    $("#drop_height_in").value = b.drop_height_in;
    $("#winglet_h_m").value = (+(b.winglet_h_m || 0)).toFixed(4);
    $("#tail_vol_h").value = (+(b.tail_vol_h || 0)).toFixed(3);
    $("#wing_prop_layout").value = b.wing_prop_layout || "mid";
    packNote();
    await evaluate();
    document.querySelector('[data-t="score"]').click();
  }
}

let catTimer = null;
async function loadCatalogue(){
  const q = encodeURIComponent($("#catq").value || "");
  const j = await (await fetch("/api/catalogue?q=" + q)).json();
  $("#cat").textContent = j.text;
}

(async function init(){
  buildTabs();
  META = await (await fetch("/api/meta")).json();
  buildForm(META);
  $("#go").onclick = evaluate;
  $("#search").onclick = startSearch;
  $("#imperial").onchange = evaluate;
  $("#catq").oninput = () => {
    clearTimeout(catTimer); catTimer = setTimeout(loadCatalogue, 200); };
  document.querySelectorAll("#left input, #left select")
    .forEach(el => el.addEventListener("change", evaluate));
  await evaluate();
  await loadCatalogue();
  $("#rules").textContent = (await (await fetch("/api/rules")).json()).text;
})();
</script>
</body>
</html>
"""
