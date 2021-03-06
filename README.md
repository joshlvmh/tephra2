<!--tephra2
===================
03-12-2014
Initial upload of the tephra2 probability scripts (perl scripts) which use the tephra2 model of ash dispersion to create survivor curves for specific locations around a volcano.

07-09-2013
Initial upload to GitHub. tephra2 is a tephra dispersion model based on an analytical (closed-form) solution of the advection-diffusion equation. This initial version uses a plume model based on the beta function and incorporates a 2-D wind model that varies in speed and direction with changes in elevation. 
-->
# tephra2

### Proposed changes:  

 - Compiler options  
 - Vectorisation key in tephra2\_calc.c  
 - Lookup table for exp() functions  
 - -xBroadwell / hardware specific flags - reduced portability, need to discuss.

### Vectorisation

```
reference pt->calculated_phi[i] has unaligned access  [  tephra2_calc.c(125,2)  ]  
```

/forward_src/new_tephra2.c:399:

```
 pt = (POINT *)GC_MALLOC((size_t)num_pts * sizeof(POINT));
```

/include/leak_detector.h:
```
#define GC_DEBUG
#include "gc.h"
#define malloc(n) GC_MALLOC(n)
#define calloc(m,n) GC_MALLOC((m)*(n))
#define free(p) GC_FREE(p)
#define realloc(p,n) GC_REALLOC((p),(n))
#define CHECK_LEAKS() GC_gcollect()
```

### Implemented changes:  

1. -Ofast / -ffast-math compiler flags 
/forward_src/makefile


### Runtimes: (for 2 runs)  

Initital:  00:54:51  
1.  00:27:39  



