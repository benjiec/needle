import { getFastaSequenceById, joinGenomeTaxonomy } from "./utils.js";
import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

const statusEl = document.getElementById('status');
const idEl = document.getElementById('identifiers');
const titleEl = document.getElementById('title');
const subtitleEl = document.getElementById('subtitle');
const consensusEl = document.getElementById('consensus');
const muscleEl = document.getElementById('muscle-msa');
const hmmalignEl = document.getElementById('hmmalign-msa');
const proteinsContainer = document.getElementById('proteins');

const queryString = window.location.search;
const urlParams = new URLSearchParams(queryString);

const moduleId = urlParams.get("module").toUpperCase();
const koId = urlParams.get("ko").toUpperCase();

function setStatus(message, type) {
  statusEl.textContent = message;
  statusEl.classList.toggle('error', type === 'error');
  statusEl.classList.toggle('hidden', type === 'hidden');
}

function setMsaImages() {
  muscleEl.src = moduleId.toLowerCase()+"/"+koId+"-muscle.png";
  hmmalignEl.src = moduleId.toLowerCase()+"/"+koId+"-hmmalign.png";
}

async function loadOrthologInfo() {
  const module_rows = await d3.tsv('../data/modules.tsv');
  const module_list = d3.index(module_rows, d => d['Module ID']);
  const ko_rows = await d3.tsv('../data/ko.tsv');
  const ko_list = d3.index(ko_rows, d => d['Ortholog ID']);
  const module = module_list.get(moduleId);
  const ko = ko_list.get(koId);

  var module_name;
  var ortholog_name;

  if (!module) { module_name = "Unknown Module"; }
  else { module_name = module["Module Name"]; }

  if (!ko) { ortholog_name = "Unknown Ortholog"; }
  else { ortholog_name = ko["Ortholog Name"]; }

  idEl.innerHTML = "<a href='./'>Modules</a>: "+moduleId + " > " + koId;
  titleEl.textContent = koId + ": " + ortholog_name;  
  subtitleEl.textContent = "Module: "+module_name;
}

function renderProteins(rowsIterable) {
  const data = Array.from(rowsIterable);
  const grid = new gridjs.Grid({
    columns: [
      { id: 'protein_hit_id', name: 'Hit' },
      { id: 'target_genome_accession', name: 'Genome' },
      { id: 'Organism', name: 'Organism' },
      { id: 'Genus', name: 'Genus' },
      { id: 'Family', name: 'Family' },
      { id: 'target_accession', name: 'Contig' },
      { id: 'hmmsearch_score', name: 'HMM Score' },
      { id: 'hmmsearch_evalue', name: 'HMM e-value' },
      { id: 'collated_protein_sequence', name: 'Sequence' }
    ],
    data,
    search: true,
    sort: true,
    pagination: { enabled: true, limit: 200 },
  });
  grid.render(proteinsContainer);
}

async function loadProteins() {
  const timestamp_param = "?t="+Date.now();
  const module_protein_rows = await d3.tsv('../data/'+moduleId.toLowerCase()+'_results/proteins.tsv'+timestamp_param);
  const ko_protein_rows = module_protein_rows.filter(d => d['query_accession'] === koId);
  await joinGenomeTaxonomy(ko_protein_rows, "target_genome_accession");
  renderProteins(ko_protein_rows);
}

async function loadConsensus() {
  const query = await d3.text('../data/'+moduleId.toLowerCase()+'_query.faa');
  const faa = getFastaSequenceById(query, koId);
  consensusEl.textContent = ">"+koId+"\n"+faa;
}

setStatus('Loading…');
setMsaImages();
loadOrthologInfo();
await loadConsensus();
await loadProteins();
setStatus('', 'hidden');

