import numpy as np
import tskit


def make_allele_count_list(ts: tskit.TreeSequence):
    parent = [tskit.NULL] * (ts.num_nodes + 1)
    num_samples_below = [0] * ts.num_nodes
    num_samples_with_derived_state = [0] * ts.num_nodes
    for s in ts.samples():
        num_samples_below[s] = 1
    mutation_node = ts.tables.mutations.node
    current_site_index = 0
    current_mutation_index = 0
    allele_count_list = []
    for diffs in ts.edge_diffs():
        for o in diffs.edges_out:
            raise NotImplementedError()
        right = diffs.interval.right
        for i in diffs.edges_in:
            parent[i.child] = i.parent
            num_samples_below[i.parent] += num_samples_below[i.child]

        # Advance sites to current tree
        while (
            current_site_index < ts.num_sites
            and ts.site(current_site_index).position < diffs.interval.left
        ):
            current_site_index += 1

        while (
            current_mutation_index < ts.num_mutations
            and current_site_index < ts.num_sites
            and ts.mutation(current_mutation_index).site != current_site_index
        ):
            current_mutation_index += 1

        while (
            current_site_index < ts.num_sites
            and ts.site(current_site_index).position < diffs.interval.right
        ):
            while (
                current_mutation_index < ts.num_mutations
                and ts.mutation(current_mutation_index).site != current_site_index
            ):
                current_mutation_index += 1
            first_mut_in_range = current_mutation_index
            last_mut_in_range = current_mutation_index
            while (
                last_mut_in_range < ts.num_mutations
                and ts.mutation(last_mut_in_range).site == current_site_index
            ):
                last_mut_in_range += 1
            assert (
                last_mut_in_range - first_mut_in_range >= 1
            ), f"{last_mut_in_range} - {first_mut_in_range} = {last_mut_in_range - first_mut_in_range}"
            mut_at_site = 0
            num_muts_at_site = last_mut_in_range - first_mut_in_range
            allele_counts = [ts.num_samples]
            alleles_at_site = [ts.site(current_site_index).ancestral_state]
            while mut_at_site < num_muts_at_site:
                mut_index = last_mut_in_range - mut_at_site - 1
                node = ts.mutation(mut_index).node
                temp = mut_at_site
                mindex = last_mut_in_range - temp - 1
                raise NotImplementedError("geet rid of the next if statement")
                if (
                    ts.mutation(mindex).derived_state
                    != ts.site(current_site_index).ancestral_state
                ):
                    num_samples_with_derived_state[parent[node]] += 1
                    nd = num_samples_below[node] - num_samples_with_derived_state[node]
                    if nd > 0:
                        try:
                            mut_allele = alleles_at_site.index(
                                ts.mutation(mut_index).derived_state
                            )
                        except ValueError:
                            mut_allele = None
                        if mut_allele is None:
                            alleles_at_site.append(ts.mutation(mut_index).derived_state)
                            mut_allele = len(alleles_at_site) - 1
                            assert mut_allele == alleles_at_site.index(
                                ts.mutation(mut_index).derived_state
                            )
                            print(
                                f"adding allele {ts.mutation(mut_index).derived_state} at {mut_allele}"
                            )
                            allele_counts.append(0)
                        allele_counts[mut_allele] += nd
                        allele_counts[0] -= nd
                temp += 1
                while (
                    temp < num_muts_at_site
                    and ts.mutation(last_mut_in_range - temp - 1).node == node
                ):
                    mindex = last_mut_in_range - temp - 1
                    temp += 1
                mut_at_site = temp
            print(f"alleles = {alleles_at_site}")
            allele_count_list.append(allele_counts)

            current_site_index += 1
            current_mutation_index = last_mut_in_range

    assert current_site_index == ts.num_sites
    assert (
        current_mutation_index == ts.num_mutations
    ), f"{current_mutation_index} != {ts.num_mutations}"
    return allele_count_list


def total_diversity(allele_count_list, n):
    div = 0.0
    for ac in allele_count_list:
        homozygosity = 0.0
        for i in ac:
            homozygosity += i * (i - 1) / (n * (n - 1))
        div += 1.0 - homozygosity
    return div


def test_case_0():
    tables = tskit.TableCollection(10.0)

    n0 = tables.nodes.add_row(0, time=2)
    n1 = tables.nodes.add_row(0, time=1)
    n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0)
    n3 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0)
    n4 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0)

    tables.edges.add_row(0.0, tables.sequence_length, n1, n2)
    tables.edges.add_row(0.0, tables.sequence_length, n1, n3)
    tables.edges.add_row(0.0, tables.sequence_length, n0, n1)
    tables.edges.add_row(0.0, tables.sequence_length, n0, n4)

    site = tables.sites.add_row(5, ancestral_state="A")
    tables.mutations.add_row(site, node=n1, time=1.0, derived_state="G")
    tables.mutations.add_row(site, node=n3, time=0.1, derived_state="C")
    tables.mutations.add_row(site, node=n1, time=1.1, derived_state="T")

    tables.sort()

    ts = tables.tree_sequence()

    print(ts.draw_text())

    allele_count_list = make_allele_count_list(ts)
    assert len(allele_count_list) == 1
    assert allele_count_list[0] == [1, 1, 1]

    div = total_diversity(allele_count_list, ts.num_samples)
    assert div == 1.0


# FIXME: the pathogical case taht we do not handle properly
# is branches to the same derived_state on different branches
def test_case_1():
    tables = tskit.TableCollection(10.0)

    n0 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n1 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n3 = tables.nodes.add_row(0, time=1.0)
    n4 = tables.nodes.add_row(0, time=2.0)

    tables.edges.add_row(0, 10, 4, 3)
    tables.edges.add_row(0, 10, 4, 2)
    tables.edges.add_row(0, 10, 3, 1)
    tables.edges.add_row(0, 10, 3, 0)

    s0 = tables.sites.add_row(5.0, ancestral_state="T")

    m0 = tables.mutations.add_row(s0, node=n0, time=0.5, derived_state="G")
    m1 = tables.mutations.add_row(s0, node=n1, time=0.5, derived_state="G")
    m2 = tables.mutations.add_row(s0, node=n2, time=0.5, derived_state="G")
    m3 = tables.mutations.add_row(s0, node=n3, time=1.5, derived_state="A")

    tables.sort()
    print(tables.mutations)
    ts = tables.tree_sequence()
    print(ts.draw_text())
    allele_counts = make_allele_count_list(ts)
    assert len(allele_counts) == 1
    assert allele_counts[0] == [0, 3]


def test_case_2():
    tables = tskit.TableCollection(10.0)

    n0 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n1 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n3 = tables.nodes.add_row(0, time=1.0)
    n4 = tables.nodes.add_row(0, time=2.0)

    tables.edges.add_row(0, 10, 4, 3)
    tables.edges.add_row(0, 10, 4, 2)
    tables.edges.add_row(0, 10, 3, 1)
    tables.edges.add_row(0, 10, 3, 0)

    s0 = tables.sites.add_row(5.0, ancestral_state="A")

    m0 = tables.mutations.add_row(s0, node=n0, time=0.5, derived_state="C")
    m1 = tables.mutations.add_row(s0, node=n1, time=0.5, derived_state="C")
    m2 = tables.mutations.add_row(s0, node=n2, time=1.5, derived_state="C")
    m3 = tables.mutations.add_row(s0, node=n3, time=1.5, derived_state="G")

    tables.sort()
    print(tables.mutations)
    ts = tables.tree_sequence()
    print(ts.draw_text())
    for i in ts.haplotypes():
        print(i)
    allele_counts = make_allele_count_list(ts)
    assert len(allele_counts) == 1
    assert allele_counts[0] == [0, 3]


def test_case_3():
    tables = tskit.TableCollection(10.0)

    n0 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n1 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n3 = tables.nodes.add_row(0, time=1.0)
    n4 = tables.nodes.add_row(0, time=2.0)

    tables.edges.add_row(0, 10, 4, 3)
    tables.edges.add_row(0, 10, 4, 2)
    tables.edges.add_row(0, 10, 3, 1)
    tables.edges.add_row(0, 10, 3, 0)

    s0 = tables.sites.add_row(5.0, ancestral_state="A")

    m0 = tables.mutations.add_row(s0, node=n0, time=0.5, derived_state="C")
    m2 = tables.mutations.add_row(s0, node=n2, time=1.5, derived_state="C")
    m3 = tables.mutations.add_row(s0, node=n3, time=1.5, derived_state="G")

    tables.sort()
    print(tables.mutations)
    ts = tables.tree_sequence()
    print(ts.draw_text())
    for i in ts.haplotypes():
        print(i)
    allele_counts = make_allele_count_list(ts)
    assert len(allele_counts) == 1
    assert allele_counts[0] == [0, 2, 1]


def test_case_4():
    tables = tskit.TableCollection(10.0)

    n0 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n1 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0.0)
    n3 = tables.nodes.add_row(0, time=1.0)
    n4 = tables.nodes.add_row(0, time=2.0)

    tables.edges.add_row(0, 10, 4, 3)
    tables.edges.add_row(0, 10, 4, 2)
    tables.edges.add_row(0, 10, 3, 1)
    tables.edges.add_row(0, 10, 3, 0)

    s0 = tables.sites.add_row(5.0, ancestral_state="A")

    m0 = tables.mutations.add_row(s0, node=n0, time=0.5, derived_state="C")
    m2 = tables.mutations.add_row(s0, node=n2, time=1.5, derived_state="C")
    m3 = tables.mutations.add_row(s0, node=n3, time=1.5, derived_state="A")

    tables.sort()
    print(tables.mutations)
    ts = tables.tree_sequence()
    print(ts.draw_text())
    for i in ts.haplotypes():
        print(i)
    allele_counts = make_allele_count_list(ts)
    assert len(allele_counts) == 1
    assert allele_counts[0] == [1, 2]
    # raise NotImplementedError("need test of anc -> derived -> anc (on just 1 branch)")
