
import numpy as np

def KL01(var:int) -> int:
    return var + 5


def RUN_SCRIPT() -> np.array:
    """Calculate krp_tc (1/d).

    Args:
        TwaterC: Water temperature (C)
        krp: Algal respiration rate at 20 degree (1/d)
    """
    to_return=0
    for i in range(1,10) :

        to_return = KL01(to_return)
        print(to_return)
        return to_return


if __name__ == "__main__":
    RUN_SCRIPT()