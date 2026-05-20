"""Custom exceptions for REFOLD."""


class REFOLDError(Exception):
    """Base exception for all REFOLD errors."""


class StructureNotFoundError(REFOLDError):
    """Raised when AlphaFold structure cannot be retrieved."""

    def __init__(self, uniprot_id: str) -> None:
        super().__init__(f"No AlphaFold structure found for UniProt ID: {uniprot_id}")
        self.uniprot_id = uniprot_id


class InvalidMutationError(REFOLDError):
    """Raised when a mutation specification is invalid."""

    def __init__(self, message: str, position: int | None = None) -> None:
        super().__init__(message)
        self.position = position


class PocketDetectionError(REFOLDError):
    """Raised when fpocket fails or produces no pockets."""


class MoleculeGenerationError(REFOLDError):
    """Raised when diffusion model fails to generate valid molecules."""


class ModelNotFoundError(REFOLDError):
    """Raised when a required model checkpoint is missing."""

    def __init__(self, checkpoint_path: str) -> None:
        super().__init__(f"Model checkpoint not found: {checkpoint_path}")
        self.checkpoint_path = checkpoint_path


class DataDownloadError(REFOLDError):
    """Raised when a required data file cannot be downloaded."""


class FilteringError(REFOLDError):
    """Raised when molecule filtering encounters an unrecoverable error."""


class SplitLeakageError(REFOLDError):
    """Raised when train/val/test splits have sequence identity leakage."""
