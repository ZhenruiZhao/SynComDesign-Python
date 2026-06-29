# Medium and Cross-Feeding Rules

The Python implementation follows the fixed MATLAB community-medium rule.

## Reaction Classes

- `external_medium_exchange`: external shared exchange reactions such as `R_EX_no3_u`.
- `strain_shared_interface`: strain exchange reactions connected to shared metabolites.
- `internal_transport`: transport reactions inside a strain model.
- `metabolic_reaction`: ordinary metabolic reactions.
- `unknown`: reactions that do not match the known classes.

Only `external_medium_exchange` can be modified by medium handling.

## Medium Application

1. Find all external shared exchange reactions.
2. Set unlisted external shared uptake lower bounds to `0`.
3. Apply listed medium lower/upper bounds to mapped external shared exchanges.
4. Under anaerobic conditions, set only external oxygen uptake to `0`.
5. Under microaerobic conditions, restrict only external oxygen uptake.
6. Do not modify strain-interface reactions.
7. Do not modify internal transport reactions.
8. Do not modify metabolic reactions.
9. Write missing medium mappings to `medium_mapping_warnings.tsv`.

## Cross-Feeding

The medium closes only outside input. A metabolite secreted by strain A into the shared pool can still be consumed by strain B through its strain-interface reaction. This is allowed because the shared metabolite participates in mass balance.

