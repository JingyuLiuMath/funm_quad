# funm

Please **do not** modify any original function. One should copy it and create a new function, then modify it.

## functions

- `sarnoldi`: sketched Arnoldi
- `sarnoldi_last_orth`: sketched Arnoldi with a  rank-1 update such that the last vector $v_{m + 1} \perp V_{m}$.
- `tarnoldi_last_orth`: truncated Arnoldi with a  rank-1 update such that the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_fom_last_orth_tarnoldi`: fom method using t-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_fom_last_orth_sarnoldi`: fom method using s-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $v_{m + 1} \perp V_{m}$.
- `funm_quad_sfom_last_sorth_tarnoldi`: sfom method using t-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $S v_{m + 1} \perp S V_{m}$.
- `funm_quad_sfom_last_sorth_sarnoldi`: sfom method using s-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $S v_{m + 1} \perp S V_{m}$.
- `funm_quad_ada_fom_last_orth_tarnoldi`: adaptive fom method using t-Arnoldi to obtain a Arnoldi-like decomposition but the last vector $S v_{m + 1} \perp S V_{m}$.


## demos

- `demo_xxx_fom_last_orth`:  using fom with last col orth to compute xxx.
- `demo_xxx_sfom_last_sorth`:  using fsom with last col sorth to compute xxx
- `demo_xxx_all_methods`: test on all methods on xxx