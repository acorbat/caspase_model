from pysb import *
from pysb.core import *
from pysb.macros import *
from pysb.macros import _macro_rule, _verify_sites, _verify_sites_complex


def cleave_complex(enzyme, e_site, substrate, s_site, product, klist, m1=None, m2=None):
    """Generate the two-step cleaving reaction
    E + S:S2 | E:S:S2 >> E + S + S2, while allowing complexes to serve as
    enzyme, substrate and/or product.

    E:S1 + S:S2 | E:S1:S:S2 >> E:S1 + S + S2

    Parameters
    ----------
    enzyme, substrate, product : Monomer, MonomerPattern, or ComplexPattern
        Monomers or complexes participating in the binding reaction.

    e_site, s_site : string
        The names of the sites on `enzyme` and `substrate` (respectively) where
        they bind each other to form the E:S complex.

    klist : list of 3 Parameters or list of 3 numbers
        Forward, reverse and cleaving rate constants (in that order). If
        Parameters are passed, they will be used directly in the generated
        Rules. If numbers are passed, Parameters will be created with
        automatically generated names based on the names and states of enzyme,
        substrate and product and these parameters will be included at the end
        of the returned component list.

    m1, m2 : Monomer or MonomerPattern
        If enzyme or substrate binding site is present in multiple monomers
        within a complex, the specific monomer desired for binding must be
        specified.

    Notes
    -----
    Binding site between enzyme and substrate is using edge 50, so this
    edge must be free to use.

    Returns
    -------
    components : ComponentSet
        The generated components. Contains the bidirectional binding Rule
        and optionally three Parameters if klist was given as numbers.
    """
    if isinstance(m1, Monomer):
        m1 = m1()
    if isinstance(m2, Monomer):
        m2 = m2()

    def build_complex(s1, site1, m1):
        _verify_sites_complex(s1, site1)
        # Retrieve a dictionary specifying the MonomerPattern within the
        # complex that contains the given binding site.
        specsitesdict = _verify_sites_complex(s1, site1)
        s1complexpatub, s1complexpatb = check_sites_comp_build(
            s1, site1, m1, specsitesdict
        )
        return s1complexpatb, s1complexpatub

    def check_sites_comp_build(s1, site1, m1, specsitesdict):
        # Return error if binding site exists on multiple monomers and a
        # monomer for binding (m1) hasn't been specified.
        if len(specsitesdict) > 1 and m1 is None:
            raise ValueError(
                f"Binding site {site1} present in more than one monomer in complex {s1}."  # noqa: E501
                "Specify variable m1, the monomer used for binding within the complex."
            )
        if not s1.is_concrete:
            raise ValueError("Complex '%s' must be concrete." % (s1))
            # If the given binding site is only present in one monomer in the complex:
        if m1 is None:
            # Build up ComplexPattern for use in rule
            # (with state of given binding site specified).
            s1complexpatub = list(specsitesdict.keys())[0]({site1: None})
            s1complexpatb = list(specsitesdict.keys())[0]({site1: 50})
            for monomer in s1.monomer_patterns:
                if monomer not in specsitesdict.keys():
                    s1complexpatub %= monomer
                    s1complexpatb %= monomer

        # If the binding site is present on more than one monomer in the
        # complex, the monomer must be specified by the user.  Use specified m1
        # to build ComplexPattern.
        else:
            # Make sure binding states of MonomerPattern m1 match those of the
            # monomer within the ComplexPattern s1 (ComplexPattern monomer
            # takes precedence if not).
            i = 0
            identical_monomers = []
            for mon in s1.monomer_patterns:
                # Only change the binding site for the first monomer that
                # matches.  Keep any others unchanged to add to final complex
                # that is returned.
                if mon.monomer.name == m1.monomer.name:
                    i += 1
                    if i == 1:
                        s1complexpatub = m1({site1: None})
                        s1complexpatb = m1({site1: 50})
                    else:
                        identical_monomers.append(mon)
            # Build up ComplexPattern for use in rule (with state of given
            # binding site  on m1 specified).
            for mon in s1.monomer_patterns:
                if mon.monomer.name != m1.monomer.name:
                    s1complexpatub %= mon
                    s1complexpatb %= mon
            if identical_monomers:
                for i in range(len(identical_monomers)):
                    s1complexpatub %= identical_monomers[i]
                    s1complexpatb %= identical_monomers[i]

        return s1complexpatub, s1complexpatb

    # If no complexes exist in the reaction, revert to catalyze().
    if (isinstance(enzyme, MonomerPattern) or isinstance(enzyme, Monomer)) and (
        isinstance(substrate, MonomerPattern) or isinstance(substrate, Monomer)
    ):
        _verify_sites(enzyme, e_site)
        _verify_sites(substrate, s_site)
        return catalyze(
            enzyme,
            e_site,
            substrate,
            s_site,
            product,
            klist,
        )

    # Build E:S
    if isinstance(enzyme, ComplexPattern):
        enzymepatb, enzyme_free = build_complex(enzyme, e_site, m1)
    else:
        enzymepatb, enzyme_free = enzyme({e_site: 50}), enzyme({e_site: None})

    if isinstance(substrate, ComplexPattern):
        substratepatb, substratepatub = build_complex(substrate, s_site, m2)
    else:
        substratepatb = substrate({s_site: 50})

    es_complex = enzymepatb % substratepatb

    if isinstance(product, ReactionPattern):
        final_product = product
        final_product += enzyme_free
    else:
        final_product = enzyme_free + product

    # Use bind complex to binding rule.

    components = bind_complex(enzyme, e_site, substrate, s_site, klist[0:2], m1, m2)
    components |= _macro_rule("cleave", es_complex >> final_product, [klist[2]], ["kc"])
    return components


def cleave_dimer(enzyme, e_site, substrate, s_site, c_site, klist):
    """Simplified macro for dimer cleaving where the reaction is
    E + M:M | E:M:M > E + M + M, and M:M is the substrate. It will generate the
    product from cleaving the substrate into it's monomer components. If both
    monomers are the same, then it will assume the first one binds to the
    enzyme.

    Parameters
    ----------
    enzyme, substrate, product : Monomer, MonomerPattern, or ComplexPattern
        Monomers or complexes participating in the binding reaction.

    e_site, s_site , c_site : string
        The names of the sites on `enzyme` and `substrate` (respectively) where
        they bind each other to form the E:S complex. The name of the site that
        is cleaved by the enzyme.

    klist : list of 3 Parameters or list of 3 numbers
        Forward, reverse and cleaving rate constants (in that order). If
        Parameters are passed, they will be used directly in the generated
        Rules. If numbers are passed, Parameters will be created with
        automatically generated names based on the names and states of enzyme,
        substrate and product and these parameters will be included at the end
        of the returned component list.

    Notes
    -----
    Binding site between enzyme and substrate is using edge 50, so this
    edge must be free to use.
    If substrate is composed of two identical monomers, the first one is
    assumed to bind the enzyme.

    Returns
    -------
    components : ComponentSet
        The generated components. Contains the bidirectional binding Rule
        and optionally three Parameters if klist was given as numbers.
    """
    # Decompose substrate to generate products
    if not isinstance(substrate, ComplexPattern):
        raise ValueError("Substrate is not a ComplexPattern")

    # If substrate is homodimer, choose the first one as binding site to enzyme
    substrate_composition = substrate.monomer_patterns
    if substrate_composition[0].monomer is substrate_composition[1].monomer:
        m2 = substrate_composition[0]
    else:
        m2 = None

    # Generate product ReactionPattern
    product = substrate_composition[0]({c_site: None})
    product += substrate_composition[1]({c_site: None})

    return cleave_complex(enzyme, e_site, substrate, s_site, product, klist, m2=m2)
