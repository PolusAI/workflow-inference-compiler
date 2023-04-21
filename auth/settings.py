import logging
from pathlib import Path
from typing import Any, Dict, Optional

from pydantic import BaseSettings

# Setup logger
logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger("wic.settings")


# Load environment variables
class Settings(BaseSettings):
    """Settings for the API.
    This defines the required inputs for the environment file.
    """

    DATA_PATH: Path = Path()
    SELECTION_PATH: Path = Path()
    LOG_LEVEL: int = logging.WARNING
    ALGORITHMS: Optional[str] = None
    AUTH_BASE_URL: Optional[str] = None
    ME_ENDPOINT: Optional[str] = None
    JWKS_ENDPOINT: Optional[str] = None
    OFFLINE_USER: Optional[str] = None

    class Config:  # NOQA: D106

        env_file = ".env"


SETTINGS = Settings()
if not (isinstance(SETTINGS.DATA_PATH, Path) and SETTINGS.DATA_PATH.exists()):
    logger.warning(
        "DATA_PATH does not exist. "
        + "Make sure to create a .env file and specify a folder for serving data."
    )
if not (isinstance(SETTINGS.SELECTION_PATH, Path) and SETTINGS.SELECTION_PATH.exists()):
    logger.warning(
        "SELECTION_PATH does not exist. "
        + "Make sure to create a .env file and specify a folder for serving data."
    )

# Define static types
CLASS_TYPES = ["string"]

SUPPORTED_FORMATS = [".feather", ".arrow", ".csv", ".parquet", ".hdf5"]

STORAGE_TYPES = [".json"]

CACHED_FILES: Dict[str, Dict[str, Any]] = {}


# Reusable Constants
DEFAULT_USER: str = "default_user"
