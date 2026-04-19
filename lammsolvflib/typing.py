"""Type aliases for array arguments."""
from typing import Union

import numpy as np
from numpy.typing import NDArray

Array1D = NDArray[np.float64]
Array2D = NDArray[np.float64]
ArrayLike1D = Union[list, Array1D]
