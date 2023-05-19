#' Deep Mutational Landscape data from 28 standardised deep mutational scans
#'
#' A dataset containing combined results from 28 deep mutational scanning
#' studies. The fitness scores are transformed into the same format and
#' normalised.
#'
#' @name deep_landscape
#' @docType data
#' @source \href{https://www.embopress.org/doi/full/10.15252/msb.202110305}{Dunham and Beltrao (2021)}
#' @keywords data
#' @format A data frame with 6357 rows
#' \describe{
#'   \item{cluster}{Amino acid subtype}
#'   \item{study}{Source study}
#'   \item{gene}{Gene tested}
#'   \item{position}{Position in gene}
#'   \item{wt}{Wild type amino acid}
#'   \item{A}{Normalised fitness score when mutating to alanine}
#'   \item{C}{Normalised fitness score when mutating to cysteine}
#'   \item{D}{Normalised fitness score when mutating to aspartate}
#'   \item{E}{Normalised fitness score when mutating to glutamate}
#'   \item{F}{Normalised fitness score when mutating to phenylalanine}
#'   \item{G}{Normalised fitness score when mutating to glycine}
#'   \item{H}{Normalised fitness score when mutating to histidine}
#'   \item{I}{Normalised fitness score when mutating to isoleucine}
#'   \item{K}{Normalised fitness score when mutating to lysine}
#'   \item{L}{Normalised fitness score when mutating to leucine}
#'   \item{M}{Normalised fitness score when mutating to methionine}
#'   \item{N}{Normalised fitness score when mutating to asparagine}
#'   \item{P}{Normalised fitness score when mutating to proline}
#'   \item{Q}{Normalised fitness score when mutating to glutamine}
#'   \item{R}{Normalised fitness score when mutating to arginine}
#'   \item{S}{Normalised fitness score when mutating to serine}
#'   \item{T}{Normalised fitness score when mutating to threonine}
#'   \item{V}{Normalised fitness score when mutating to valine}
#'   \item{W}{Normalised fitness score when mutating to tryptophan}
#'   \item{Y}{Normalised fitness score when mutating to tyrosine}
#'   \item{mean_score}{Mean normalised fitness score}
#'   \item{mean_sift}{Mean SIFT4G score for the position}
#'   \item{total_energy}{Mean FoldX total ddG for mutations to this position}
#'   \item{backbone_hbond}{Mean FoldX backbone_hbond term}
#'   \item{sidechain_hbond}{Mean FoldX sidechain_hbond term}
#'   \item{van_der_waals}{Mean FoldX van_der_waals term}
#'   \item{electrostatics}{Mean FoldX electrostatics term}
#'   \item{solvation_polar}{Mean FoldX solvation_polar term}
#'   \item{solvation_hydrophobic}{Mean FoldX solvation_hydrophobic term}
#'   \item{van_der_waals_clashes}{Mean FoldX van_der_waals_clashes term}
#'   \item{entropy_sidechain}{Mean FoldX entropy_sidechain term}
#'   \item{entropy_mainchain}{Mean FoldX entropy_mainchain term}
#'   \item{cis_bond}{Mean FoldX cis_bond term}
#'   \item{torsional_clash}{Mean FoldX torsional_clash term}
#'   \item{backbone_clash}{Mean FoldX backbone_clash term}
#'   \item{helix_dipole}{Mean FoldX helix_dipole term}
#'   \item{disulfide}{Mean FoldX disulfide term}
#'   \item{electrostatic_kon}{Mean FoldX electrostatic_kon term}
#'   \item{partial_covalent_bonds}{Mean FoldX partial_covalent_bonds term}
#'   \item{energy_ionisation}{Mean FoldX energy_ionisation term}
#'   \item{phi}{Phi backbone angle}
#'   \item{psi}{Psi backbone angle}
#'   \item{all_atom_rel}{Relative surface accessibility}
#'   \item{hydrophobicity}{Residue hydrophobicity}
#'   \item{PC1}{PC1 of the principal component analysis of A-Y}
#'   \item{PC2}{PC2 of the principal component analysis of A-Y}
#'   \item{PC3}{PC3 of the principal component analysis of A-Y}
#'   \item{PC4}{PC4 of the principal component analysis of A-Y}
#'   \item{PC5}{PC5 of the principal component analysis of A-Y}
#'   \item{PC6}{PC6 of the principal component analysis of A-Y}
#'   \item{PC7}{PC7 of the principal component analysis of A-Y}
#'   \item{PC8}{PC8 of the principal component analysis of A-Y}
#'   \item{PC9}{PC9 of the principal component analysis of A-Y}
#'   \item{PC10}{PC10 of the principal component analysis of A-Y}
#'   \item{PC11}{PC11 of the principal component analysis of A-Y}
#'   \item{PC12}{PC12 of the principal component analysis of A-Y}
#'   \item{PC13}{PC13 of the principal component analysis of A-Y}
#'   \item{PC14}{PC14 of the principal component analysis of A-Y}
#'   \item{PC15}{PC15 of the principal component analysis of A-Y}
#'   \item{PC16}{PC16 of the principal component analysis of A-Y}
#'   \item{PC17}{PC17 of the principal component analysis of A-Y}
#'   \item{PC18}{PC18 of the principal component analysis of A-Y}
#'   \item{PC19}{PC19 of the principal component analysis of A-Y}
#'   \item{PC20}{PC20 of the principal component analysis of A-Y}
#'   \item{tSNE1}{First tSNE dimension of dimensionality reduction of A-Y}
#'   \item{tSNE2}{Second tSNE dimension of dimensionality reduction of A-Y}
#'   \item{umap1}{First UMAP dimension of dimensionality reduction of A-Y}
#'   \item{umap2}{Second UMAP dimension of dimensionality reduction of A-Y}
#' }
"deep_landscape"

#' Amino acid positional subtypes from the deep mutational landscape
#'
#' A dataset containing groups, descriptions and notes about each amino acid subtype identified in the deep mutational
#' landscape. Groups were determined by ER profile correlation and descriptions/notes from manual assessment of subtype
#' position properties (ER, FoldX predicitons, conservation, etc.).
#'
#' @name subtypes
#' @docType data
#' @source \href{https://www.embopress.org/doi/full/10.15252/msb.202110305}{Dunham and Beltrao (2021)}
#' @keywords data
#' @format A data frame with 119 rows and 5 columns
#' \describe{
#'   \item{wt}{Wild type amino acid}
#'   \item{cluster}{Amino acid subtype}
#'   \item{prop}{Proportion of amino acid positions in the original dataset that are classified into this cluster}
#'   \item{group}{Subtype group}
#'   \item{description}{Description of the subtype}
#'   \item{notes}{Additional notes about positions of this subtype}
#'   \item{A}{Normalised fitness score when mutating to alanine}
#'   \item{C}{Normalised fitness score when mutating to cysteine}
#'   \item{D}{Normalised fitness score when mutating to aspartate}
#'   \item{E}{Normalised fitness score when mutating to glutamate}
#'   \item{F}{Normalised fitness score when mutating to phenylalanine}
#'   \item{G}{Normalised fitness score when mutating to glycine}
#'   \item{H}{Normalised fitness score when mutating to histidine}
#'   \item{I}{Normalised fitness score when mutating to isoleucine}
#'   \item{K}{Normalised fitness score when mutating to lysine}
#'   \item{L}{Normalised fitness score when mutating to leucine}
#'   \item{M}{Normalised fitness score when mutating to methionine}
#'   \item{N}{Normalised fitness score when mutating to asparagine}
#'   \item{P}{Normalised fitness score when mutating to proline}
#'   \item{Q}{Normalised fitness score when mutating to glutamine}
#'   \item{R}{Normalised fitness score when mutating to arginine}
#'   \item{S}{Normalised fitness score when mutating to serine}
#'   \item{T}{Normalised fitness score when mutating to threonine}
#'   \item{V}{Normalised fitness score when mutating to valine}
#'   \item{W}{Normalised fitness score when mutating to tryptophan}
#'   \item{Y}{Normalised fitness score when mutating to tyrosine}
#'   \item{mean_score}{Mean normalised ER fitness score}
#'   \item{mean_sift}{Mean SIFT4G score for positions of this cluster}
#'   \item{total_energy}{Mean FoldX total ddG for mutations to position of this cluster}
#'   \item{backbone_hbond}{Mean FoldX backbone_hbond term}
#'   \item{sidechain_hbond}{Mean FoldX sidechain_hbond term}
#'   \item{van_der_waals}{Mean FoldX van_der_waals term}
#'   \item{electrostatics}{Mean FoldX electrostatics term}
#'   \item{solvation_polar}{Mean FoldX solvation_polar term}
#'   \item{solvation_hydrophobic}{Mean FoldX solvation_hydrophobic term}
#'   \item{van_der_waals_clashes}{Mean FoldX van_der_waals_clashes term}
#'   \item{entropy_sidechain}{Mean FoldX entropy_sidechain term}
#'   \item{entropy_mainchain}{Mean FoldX entropy_mainchain term}
#'   \item{cis_bond}{Mean FoldX cis_bond term}
#'   \item{torsional_clash}{Mean FoldX torsional_clash term}
#'   \item{backbone_clash}{Mean FoldX backbone_clash term}
#'   \item{helix_dipole}{Mean FoldX helix_dipole term}
#'   \item{disulfide}{Mean FoldX disulfide term}
#'   \item{electrostatic_kon}{Mean FoldX electrostatic_kon term}
#'   \item{partial_covalent_bonds}{Mean FoldX partial_covalent_bonds term}
#'   \item{energy_ionisation}{Mean FoldX energy_ionisation term}
#'   \item{all_atom_rel}{Mean relative surface accessibility}
#' }
"subtypes"

#' Example deep_mutational_scan objects
#'
#' A list of \code{\link{deep_mutational_scan}} objects. This dataset is primarily used to illustrate use of
#' the package functions. The studies were sourced from \url{https://www.mavedb.org}{MaveDB}.
#'
#' @name deep_scans
#' @docType data
#' @keywords data
#' @format A list of \code{\link{deep_mutational_scan}} objects
#' \describe{
#'   \item{hsp90}{HSP90 scan by Hietpas et al. (2011, doi: 10.1073/pnas.1016024108, MaveDB: urn:mavedb:00000011-a-1)}
#'   \item{gpa}{GpA scan by Elazar et al. (2016, doi: 10.7554/eLife.12125, MaveDB:  urn:mavedb:00000051-c-1)}
#'   \item{p53}{p53 scan by Kotler et al. (2018, doi: 10.1016/j.molcel.2018.06.012 , MaveDB: urn:mavedb:00000059-c-1)}
#' }
"deep_scans"
