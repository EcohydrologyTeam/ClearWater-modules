from abc import ABC, abstractmethod
from custom_types import ArrayLike
from pathlib import Path
from datetime import datetime, timedelta


class DataProvider(ABC):
    """
    Base class for providers.
    """

    @abstractmethod
    def __init__(self, **kwargs) -> None:
        raise NotImplementedError

    @abstractmethod
    def read(self) -> ArrayLike:
        raise NotImplementedError

    @abstractmethod
    def write_to_store(
        self,
        store_path: Path,
        start_time: datetime,
        end_time: datetime,
        time_step: timedelta,
        variable_name: str,
        field_name: str | None = None,
    ) -> None:
        raise NotImplementedError
