import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import graph_tool.all as gt



def adata_to_graph(
    adata : ad.AnnData, 
    obs_list : str | list[str] | None = None,
    weight_name : str = "weight",
    categorical_as_codes : bool = False,
    add_colors : bool = False
):
    """
    Convert AnnData connectivity (adata.obsp['connectivities']) to a graph-tool Graph.
    Creates vertex properties from columns listed in obs_cols.

    Args:
        adata: AnnData with .obsp['connectivities'] and .obs
        obs_cols: list[str] of columns in adata.obs to convert to vertex properties.
                  If None, no vertex properties are created.
        weight_name: name for edge weight property (if adjacency has data)
        categorical_as_codes: if True, pandas.Categorical columns are stored as integer codes;
                              otherwise stored as string labels.
        add_colors: if True, adds a vertex property '{col}_color' with colors from adata.uns['{col}_colors'];
                     Only for those columns that have a colors listed in data.uns. 

    Returns:
        g: graph_tool.Graph with vertex properties added under g.vertex_properties[colname]
           and (if present) an edge property g.edge_properties[weight_name].
    """

    if isinstance(obs_list, str):
        obs_list = [obs_list]

    connectivities = adata.obsp.get('connectivities', None)
    if connectivities is None:
        raise KeyError("adata.obsp does not contain 'connectivities'")

    # ensure sparse COO
    if not isinstance(connectivities, sp.coo_matrix):
        connectivities = sp.coo_matrix(connectivities)

    directed = (connectivities != connectivities.T).nnz != 0

    n = int(adata.n_obs)
    g = gt.Graph(directed=directed)
    g.add_vertex(n)

    
    # add edge weights if adjacency contains non-1 data
    has_weights = connectivities.data is not None and np.any(connectivities.data != 1)
    if has_weights:
        eweight = g.new_edge_property("double")
        g.edge_properties[weight_name] = eweight
    else:
        # still create a weight property if you want (optional)
        eweight = None

    # Add each undirected edge once (avoid duplicates from symmetric matrix)
    rows = connectivities.row.astype(int)
    cols = connectivities.col.astype(int)
    data = connectivities.data if connectivities.data is not None else np.ones_like(rows, dtype=float)

    mask = rows < cols
    for r, c, w in zip(rows[mask], cols[mask], data[mask]):
        e = g.add_edge(int(r), int(c))
        if has_weights:
            g.edge_properties[weight_name][e] = float(w)

    # Adding the index name as a vertex property
    index_prop = g.new_vertex_property("string")
    for i, v in enumerate(g.vertices()):
        index_prop[v] = str(adata.obs.index[i])
    g.vertex_properties['index'] = index_prop
    
    # create vertex properties from obs columns
    if obs_list:
        for col in obs_list:
            if col not in adata.obs.columns:
                raise KeyError(f"{col} not found in adata.obs")
            s = adata.obs[col]
            # choose property type
            if pd.api.types.is_integer_dtype(s.dtype):
                prop = g.new_vertex_property("int", name=col)
                for i, val in enumerate(s.values):
                    prop[g.vertex(i)] = int(val) if pd.notna(val) else 0
            elif pd.api.types.is_float_dtype(s.dtype):
                prop = g.new_vertex_property("double")
                for i, val in enumerate(s.values):
                    prop[g.vertex(i)] = float(val) if pd.notna(val) else float("nan")
            elif isinstance(s.dtype, pd.CategoricalDtype):
                if categorical_as_codes:
                    prop = g.new_vertex_property("int")
                    codes = s.cat.codes.values
                    for i, val in enumerate(codes):
                        prop[g.vertex(i)] = int(val)
                    prop_int = None
                else:
                    prop = g.new_vertex_property("string")
                    for i, val in enumerate(s.astype(str).values):
                        prop[g.vertex(i)] = val
                    prop_int = g.new_vertex_property("int")
                    codes = s.cat.codes.values
                    for i, val in enumerate(codes):
                        prop_int[g.vertex(i)] = int(val)
            else:
                # fallback: store as string
                prop = g.new_vertex_property("string")
                for i, val in enumerate(s.astype(str).values):
                    prop[g.vertex(i)] = val

            # attach vertex property to graph
            g.vertex_properties[col] = prop
            if prop_int is not None:
                g.vertex_properties[col + "_code"] = prop_int

            # Adding colors if applicable
            if f"{col}_colors" in adata.uns.keys():
                color_dict = {}
                for _cat, _col in zip(adata.obs[col].cat.categories, adata.uns[f"{col}_colors"]):
                    color_dict[_cat] = _col
                colors = g.new_vertex_property("string")                
                for v in g.vertices():
                    _cat = g.vertex_properties[col][v]
                    _col = color_dict[_cat]
                    colors[v] = _col
                g.vertex_properties[f"{col}_color"] = colors

    return g


def assign_nested_cluster_property(nested_state, g=None, prop_name="nested_cluster", sep="|", add_level_props=True):
    """
    Attach nested cluster assignments from a graph-tool NestedBlockState to `g` as vertex properties.

    Args:
        nested_state: graph_tool.models.NestedBlockState (result of minimize_nested_blockmodel_dl)
        g: graph_tool.Graph instance. If None, will use nested_state.get_graph().
        prop_name: base name for the combined string property (and prefix for per-level props).
        sep: separator used in the combined string property between levels.
        add_level_props: if True, also add integer vertex properties for each nesting level
                         named f"{prop_name}_lvl{i}".
    Returns:
        dict of created property maps: {"combined": combined_prop, "levels": [lvl_prop0, ...]}
    """
    if g is None:
        g = nested_state.get_graph()

    levels = nested_state.get_levels()
    for nl, _level in enumerate(levels):
        num_clusters = _level.get_N()
        if num_clusters == 1:
            break
    n_levels = nl

    combined = g.new_vertex_property("string")
    level_props = []
    if add_level_props:
        for _ in range(n_levels):
            level_props.append(g.new_vertex_property("int"))

    for v in g.vertices():
        parts = []
        b = v
        for i, lvl in enumerate(levels):
            if i >= n_levels:
                break
            blocks = lvl.get_blocks()
            b = int(blocks[b])
            parts.append(str(b))
            if add_level_props:
                level_props[i][v] = b
        combined[v] = sep.join(parts)

    g.vertex_properties[prop_name] = combined
    if add_level_props:
        for i, p in enumerate(level_props):
            g.vertex_properties[f"{prop_name}_lvl{i}"] = p

    return {"combined": combined, "levels": level_props}


def _nested_cluster_graph(
        g : gt.Graph,
):
    """
    Cluster the graph using the Bayesian blockmodel.
    """
    # TODO (a3klein): add options for more extensive clustering via expansive search of the solution space
    state = gt.minimize_nested_blockmodel_dl(g)
    _ = assign_nested_cluster_property(state, g=g, prop_name="nested_cluster", sep="|", add_level_props=True)
    return state, g


def graph_tool_vertex_props_to_dataframe(g, props=None, index_col=None, drop_index_col=False, decode_bytes=True):
    """
    Convert graph-tool vertex properties to a pandas.DataFrame.

    Args:
        g: graph_tool.Graph
        props: list of property names to include (default: all vertex properties)
        index_col: name of a column to use as DataFrame index (default: None -> integer vertex index)
        drop_index_col: if True and index_col is set, drop that column from the DataFrame after setting index
        decode_bytes: decode bytes/bytearray values to str when encountered

    Returns:
        pandas.DataFrame with one row per vertex in graph-tool's vertex order.
    """
    if props is None:
        props = list(g.vertex_properties.keys())

    vlist = list(g.vertices())  # preserves vertex order
    idx = [int(v) for v in vlist]

    data = {}
    for name in props:
        if name not in g.vertex_properties:
            raise KeyError(f"vertex property '{name}' not found")
        vprop = g.vertex_properties[name]
        col = []
        for v in vlist:
            val = vprop[v]
            if decode_bytes and isinstance(val, (bytes, bytearray)):
                try:
                    val = val.decode()
                except Exception:
                    val = str(val)
            col.append(val)
        data[name] = col

    df = pd.DataFrame(data, index=idx)

    if index_col is not None:
        if index_col not in df.columns:
            raise KeyError(f"index_col '{index_col}' not found in selected properties")
        df.index = df[index_col].values
        if drop_index_col:
            df = df.drop(columns=[index_col])

    return df


def bayes_clust(
    adata = ad.AnnData,
    cols : str | list[str] | None = None,
):
    """
    Perform Bayesian clustering on the AnnData connectivity graph.

    Args:
        adata: AnnData with .obsp['connectivities'] and .obs
        cols: list of columns in adata.obs to add as vertex properties to the graph.
              If None, no vertex properties are added.

    Returns:
        adata: updated AnnData with new .obsp['nested_cluster_graph'] and .obs columns
    """
    g = adata_to_graph(
        adata,
        obs_list=cols,
        weight_name="weight",
        categorical_as_codes=False,
        add_colors=False
    )
    state, g = _nested_cluster_graph(g)
    assign_nested_cluster_property(state, g=g, prop_name="nested_cluster", sep="|", add_level_props=True)
    df = graph_tool_vertex_props_to_dataframe(g)
    df.set_index('index', drop=True, inplace=True)
    adata.obsm['bayes_clustering'] = df
    return adata