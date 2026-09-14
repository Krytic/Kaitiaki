# Notice

This submodule was authored by Hamish Jelleyman and co-written with Claude: https://github.com/hjelleyman/stars-kippenhahn

Hamish published it under the MIT license and has explicitly allowed us to use it, modify it, and ship it as part of kaitiaki. We are in the process of more fulsomely integrating it into the codebase.

This kippenhahn diagram generator can be used by instansiating a `plotfile` object and passing `engine='new'` to the `kippenhahn_diagram()` method:

```py
import kaitiaki
import matplotlib.pyplot as plt
plot = kaitiaki.file.plot()
plot.kippenhahn_diagram(engine='new')
plt.show()
```