import numpy as np

unlog_y = lambda ax, arr: 10**arr if ax=='y' else arr
unlog_x = lambda ax, arr: 10**arr if ax=='x' else arr
log10_y = lambda ax, arr: np.log10(arr) if ax=='y' else arr
log10_x = lambda ax, arr: np.log10(arr) if ax=='x' else arr
loge_y  = lambda ax, arr: np.log(arr) if ax=='y' else arr
loge_x  = lambda ax, arr: np.log(arr) if ax=='x' else arr

unlog_both = lambda ax, arr: 10**arr