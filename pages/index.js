import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

const statusEl = document.getElementById('status');
const gridContainer = document.getElementById('grid');

function setStatus(message, type) {
  statusEl.textContent = message;
  statusEl.classList.toggle('error', type === 'error');
  statusEl.classList.toggle('hidden', type === 'hidden');
}

function renderGrid(rowsIterable) {
  const data = Array.from(rowsIterable);
  const grid = new gridjs.Grid({
    columns: [
      { id: 'Module ID', name: 'Module ID' },
      { id: 'Module Name', name: 'Module Name' },
      { id: 'Ortholog ID', name: 'Ortholog ID', formatter: (cell, row) => gridjs.html("<a href='ko.html?module="+row.cells[0].data+"&ko="+cell+"'>"+cell+"</a>") },
      { id: 'Ortholog Name', name: 'Ortholog Name' },
    ],
    data,
    search: true,
    sort: true,
    pagination: { enabled: true, limit: 200 },
  });
  grid.render(gridContainer);
}

async function load() {
  const timestamp_param = "?t="+Date.now();
  const ko_rows = await d3.tsv('../data/ko.tsv');
  const ko_list = d3.index(ko_rows, d => d['Ortholog ID']);
  const all_module_rows = await d3.tsv('../data/modules.tsv');
  const all_module_list = d3.index(all_module_rows, d => d['Module ID']);
  const all_module_ko_rows = await d3.tsv('../data/module_ko.tsv');
  const show_module_rows = await d3.tsv('./modules.tsv'+timestamp_param);

  var module_kos = [];

  show_module_rows.forEach((v) => {
    var m = all_module_list.get(v['Module ID']);
    const module_ko_rows = all_module_ko_rows.filter(d => d['Module ID'] === v['Module ID']);
    if (module_ko_rows.length === 1) {
      const module_ko_list = module_ko_rows[0]['Module Definition'];
      module_ko_list.split(",").forEach((koId) => {
        var ko = {
          'Module ID': m['Module ID'], 
          'Module Name': m['Module Name'],
          'Ortholog ID': koId,
          'Ortholog Name': ko_list.get(koId)['Ortholog Name']
        };
        module_kos.push(ko);
      });
    }
  });

  module_kos.sort((a, b) => a['Ortholog ID'].localeCompare(b['Ortholog ID']));
  renderGrid(module_kos);
}

setStatus('Loading…');
await load();
setStatus('', 'hidden');
