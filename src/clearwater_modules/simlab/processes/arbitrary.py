import xsimlab as xs
from .riverine import Riverine


@xs.process
class Arbitrary:
    """
    This is an arbitrary process that does nothing but is used to test the model.
    """

    calls = xs.variable(
        dims=(),
        intent="out",
        description="Just a number indicating home many times the process has been called",
    )
    volume = xs.foreign(
        other_process_cls=Riverine,
        var_name="volume",
        intent="in",
    )

    def initialize(self) -> None:
        self._number_of_calls = 0

    @xs.runtime(args="step")
    def run_step(self, dt):
        sum = self.volume.sum()
        self._number_of_calls += 1

    def finalize_step(self):
        self.calls = self._number_of_calls
