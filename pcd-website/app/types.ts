export interface PCDEntry {
  entry_id: string;
  investigator: string;
  institution: string;
  metadata: {
    gene: string;
    uniprot: string;
    clinvar_id: string;
    disease: string;
    mutation_precursor: string;
    mutation_mature: string;
    variant_class: string;
    mechanism: string;
  };
  pocket: {
    target_conformation: string;
    n_conformations_sampled: number;
    fpocket_druggability: number;
    volume_angstrom3: number;
    alpha_sphere_count: number | string;
    center_angstrom: [number, number, number];
    wt_baseline_druggability: number;
    exceeds_threshold: boolean;
    pocket_type: string;
    dist_mutation_to_pocket_angstrom: number;
  };
  sequence: {
    fasta_header: string;
    fasta_sequence: string;
    pocket_lining_positions_precursor?: number[];
    pocket_lining_positions?: number[];
    pocket_lining_residues: string;
    pocket_lining_1letter: string;
    sequence_slice_around_pocket: string;
  };
  eij_matrix: {
    file: string;
    shape: [number, number];
    residue_labels: string[];
    values?: number[][];
  };
  chaperone: {
    smiles: string;
    iupac_name: string;
    common_name: string;
    composite_score: number;
    mw: number;
    logp: number;
    qed: number;
    sa_score: number;
    hbd: number;
    hba: number;
    tpsa: number;
    lipinski: boolean;
    veber: boolean;
    pains: boolean;
    pocket_affinity_score: number;
    binding_mode: Record<string, string>;
  };
  assets: {
    eij_matrix_npy: string;
    eij_matrix_csv?: string;
    transient_conformation_pdb: string;
    pocket_summary_json: string;
    stage3_candidates_json: string;
    pdb_structure?: string;
  };
  status: string;
}

export interface ProteomeTargets {
  total_clinvar_pathogenic_missense: number;
  total_processed: number;
  queue_remaining: number;
}

export interface PCDAtlas {
  database: string;
  version: string;
  build_date: string;
  investigator: string;
  institution: string;
  powered_by: string;
  total_entries: number;
  proteome_targets?: ProteomeTargets;
  entries: PCDEntry[];
}
