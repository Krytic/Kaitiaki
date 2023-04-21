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
        result = 0
        for sigma in self.__sigmas:
            result += np.sqrt(2/np.pi) * vk **2 / sigma**3 * np.exp(-vk**2/(2*sigma**2))
        return result

class __Verbunt(rv_continuous):
    def __init__(self, sigmas=[75,316], weights=[0.42, 1-0.42], **kwargs):
        if 'a' not in kwargs.keys():
            kwargs['a'] = 0
        if 'b' not in kwargs.keys():
            kwargs['b'] = 2000

        super().__init__(**kwargs)

        if isinstance(sigmas, (list, np.ndarray, tuple)):
            self.__sigmas = sigmas
        else:
            raise TypeError("sigmas should be iterable")

        if len(self.__sigmas) != len(weights):
            raise ValueError("sigmas and weights must be the same length")
        elif len(weights) != 2:
            raise ValueError("Only two peaks permitted for the VIC model.")
        else:
            self.__weights = weights

        if sum(self.__weights) != 1:
            # normalise to sum to 1
            sum_of_weights = sum(self.__weights)

            self,__weights = [w/sum_of_weights for w in self.__weights]


    def _pdf(self, vk):
        return np.sqrt(2/np.pi) * vk**2 * ((self.__weights[0])/(self.__sigmas[0])**3*np.exp(-(vk**2)/(2*((self.__sigmas[0])**2)))+(self.__weights[1])/(self.__sigmas[1])**3*np.exp(-(vk**2)/(2*((self.__sigmas[1])**2))))

def __Bray2018(mej, mrem):
    mns = mrem
    if mrem > 2.5:
        mns = 1.4

    return 100*mej/mrem-170*(mns/mrem)

def __BrayRichards(mej, mrem):
    mns = mrem
    if mrem > 2.5:
        mns = 1.4

    return 115*mej/mrem+15*(mns/mrem)

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


def sample(dist_name, n_samples, **kwargs):
    if not dist_name.startswith("Bray"):
        dist = _get_dist_by_name(dist_name, **kwargs)
        return dist.rvs(size=n_samples)
    else:
        if 'mej' in kwargs and 'mrem' in kwargs:
            funcs = [cls_name for cls_name, cls_obj in inspect.getmembers(sys.modules[__name__]) if inspect.isfunction(cls_obj) and cls_name[:2] == "__"]
            dist = getattr(sys.modules[__name__], f'__{dist_name}')

            return np.array([dist(**kwargs) for _ in range(n_samples)])
        else:
            raise ValueError("If you want a Bray kick I need mej and mrem.")

def pdf(dist_name, **kwargs):
    dist = _get_dist_by_name(dist_name, **kwargs)
    return dist.pdf