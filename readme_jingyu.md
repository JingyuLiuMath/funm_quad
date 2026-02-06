# funm

Please **do not** modify any original function. One should copy it and create a new function, then modify it.

## functions

- `sketched_arnoldi`: sketched Arnoldi
- `sketched_arnoldi_last_update`: sketched Arnoldi with a  rank-1 update such that the last vector $v_{m + 1} \perp V_{m}$.
- `tarnoldi_last_update`: truncated Arnoldi with a  rank-1 update such that the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_fom_last_update_sarnoldi`: fom method using t-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_fom_last_update_sarnoldi`: fom method using s-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_sfom_last_update_sarnoldi`: sfom method using t-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_sfom_last_update_sarnoldi`: sfom method using s-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $v_{m + 1} \perp V_{m}$.

## demos

- `demo_xxx_fom_last_update`:  using fom with last col update to compute xxx.
- `demo_xxx_sfom_last_update`:  using fsom with last col update to compute xxx
- `demo_xxx_all_methods`: test on all methods on xxx