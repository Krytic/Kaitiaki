import inspect
import sys

import kaitiaki
from scipy.stats import rv_continuous
import numpy as np

class __Hobbs(rv_continuous):
    "Hobbs Distribution"
    def __init__(self, sigma=265, **kwargs):
        if 'a' not in kwargs.keys():
            kwargs['a'] = 0
        if 'b' not in kwargs.keys():
            kwargs['b'] = 2000

        super().__init__(**kwargs)
        self.__sigma = sigma

    def _pdf(self, vk_h):
        sigma = self.__sigma
        return np.sqrt(2/np.pi) * vk_h **2 / sigma**3 * np.exp(-vk_h**2/(2*sigma**2))

class __MultimodalHobbs(rv_continuous):
    def __init__(self, sigmas=265, **kwargs):
        if 'a' not in kwargs.keys():
            kwargs['a'] = 0
        if 'b' not in kwargs.keys():
            kwargs['b'] = 2000

        super().__init__(**kwargs)

        if isinstance(sigmas, (list, np.ndarray, tuple)):
            self.__sigmas = sigmas
        else:
            raise TypeError("sigmas should be iterable")


    def _pdf(self, vk):
        result = np.zeros(len(vk))
        for sigma in self.__sigmas:
            result += np.sqrt(2/np.pi) * vk **2 / sigma**3 * np.exp(-vk**2/(2*sigma**2))
        return result



def _get_dist_by_name(dist_name, **kwargs):
    internal_name = f"__{dist_name[0].upper() + dist_name[1:]}"

    # Following line gets all classes defined in this file as long
    # as they start with a __ (two underscores).
    classes = [cls_name for cls_name, cls_obj in inspect.getmembers(sys.modules[__name__]) if inspect.isclass(cls_obj) and cls_name[:2] == "__"]
    if internal_name not in classes:
        raise ValueError(f"Distribution {dist_name} is not defined in kaitiaki.")

    kwargs['name'] = dist_name
    dist = getattr(sys.modules[__name__], internal_name)(**kwargs)

    return dist


def sample(dist_name, n_samples=1, **kwargs):
    dist = _get_dist_by_name(dist_name, **kwargs)
    return dist.rvs(size=n_samples)

def pdf(dist_name, **kwargs):
    dist = _get_dist_by_name(dist_name, **kwargs)
    return dist.pdf