from abc import ABC, abstractmethod
from custom_types import ArrayLike
from pathlib import Path
from datetime import datetime, timedelta
from typing import Protocol, runtime_checkable
from variables import Variable


@runtime_checkable
class DataSource(Protocol):
    @abstractmethod
    def read(self, parameter_name: str) -> Variable:
        raise NotImplementedError


@runtime_checkable
class ChunkedDataSource(Protocol):
    @abstractmethod
    def read_chunk(
        self,
        parameter_name: str,
        start_time: datetime,
        end_time: datetime,
    ) -> Variable:
        raise NotImplementedError


@runtime_checkable
class DataStore(Protocol):
    @abstractmethod
    def write(self, data: ArrayLike, parameter_name: str) -> None:
        raise NotImplementedError


@runtime_checkable
class ChunkedDataStore(Protocol):
    @abstractmethod
    def write_chunk(
        self,
        data: ArrayLike,
        parameter_name: str,
        start_time: datetime,
        end_time: datetime,
    ) -> None:
        raise NotImplementedError


# This is what is currently implemented


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


class ChunkedDataProvider(DataProvider):
    def read_chunk(self, start_time: datetime, end_time: datetime) -> ArrayLike:
        raise NotImplementedError
