import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

/**
 * Return the sequence for a given ID from FASTA text.
 * - ID defaults to the first whitespace-delimited token in the header line.
 * - Set matchFullHeader=true to match the entire header (without the leading '>').
 * - Set caseSensitive=true for case-sensitive header matching (default is case-insensitive).
 */
function getFastaSequenceById(fastaText, targetId, { matchFullHeader = false, caseSensitive = false } = {}) {
  if (!fastaText || !targetId) return null;

  const normalize = s => s.replace(/\r\n/g, '\n').replace(/\r/g, '\n');
  const lines = normalize(String(fastaText)).split('\n');

  const norm = s => (caseSensitive ? s : s.toLowerCase());
  const target = norm(targetId);

  let currentHeader = null;
  let currentSeq = [];

  const headerMatches = headerLine => {
    const headerNoGt = headerLine.slice(1).trim();
    const headerId = headerNoGt.split(/\s+/)[0]; // token before first whitespace
    const key = matchFullHeader ? headerNoGt : headerId;
    return norm(key) === target;
  };

  const flushIfMatch = () => {
    if (currentHeader && headerMatches(currentHeader)) {
      return currentSeq.join('').replace(/\s+/g, '');
    }
    return null;
  };

  for (const line of lines) {
    if (line.startsWith('>')) {
      const hit = flushIfMatch();
      if (hit !== null) return hit;
      currentHeader = line;
      currentSeq = [];
    } else if (currentHeader) {
      if (line.trim().length > 0) currentSeq.push(line.trim());
    }
  }

  return flushIfMatch();
}

/*
 * Given an array of dictionary, joins the array with genomes.tsv data, with "Genome Accession" as key.
 */
async function joinGenomeTaxonomy(dataIterable, dataKey) {
  const genome_rows = await d3.tsv('../data/genomes.tsv');
  const genomes = d3.index(genome_rows, d => d['Genome Accession']);
  dataIterable.forEach((d) => {
    const joinValue = d[dataKey];
    const genome = genomes.get(joinValue);
    if (genome) {
      Object.assign(d, genome);
    }
  })
}

export { getFastaSequenceById, joinGenomeTaxonomy };
