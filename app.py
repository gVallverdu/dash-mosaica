#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
## Documentation

This application aims to visualize structural data in order to provide a geometrical
analysis of a molecular strucutre. It also provides a way to visualize any 
atomic quantities.

The definitions of the geometrical data available by default are given below.
Some of them are available only if there is a minimum number of bonds:

* **Angular defect (degrees):** The angular defect is a measure of the discrete curvature
on a given atom. It is computed as360° minus the sum of the angles between bonds with atoms
bonded to the considered atom.
* **haddon (degrees):** This is the pyramidalization angle as defined by 
Haddon et. al. [ref]
* **improper angle (degrees):** This is the improper dihedral angle
* **dist. from. ave. plane (angstrom):** This is the distance between the
considered atom and the average plane defined by atoms bonded to it.
* **neighbors (number):** This is the number of neighbors of the atom
* **ave. neighb. dist. (angstrom):** This is the average distance between the 
considered atom and its neighbors.

### File upload

The application accepts standard xyz files. The first line of the file must
contains the number of atom, the second line is a title (not considered here) and
the following lines start by the element name followed by the cartesian coordinates.
"""

import io
import base64

import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import dash_bio

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils

__author__ = "Germain Salvato Vallverdu"
__title__ = "Structural data viewer"
__subtitle__ = "Part of the Mosaica project"

# plotly colorscales
# must be lowercase for matplotlib
COLORSCALES = ["Blues",
               "Blues_r"
               "Greens",
               "Greens_r",
               "Greys",
               "Greys_r",
               "hot",
               "hot_r",
               "jet",
               "jet_r",
               "inferno",
               "inferno_r",
               "magma",
               "magma_r",
               "plasma",
               "plasma_r",
               "rainbow",
               "rainbow_r",
               "RdBu",
               "RdBu_r",
               "Reds",
               "Reds_r",
               "viridis",
               "viridis_r",
               "cividis",
               "cividis_'",
               "YlGnBu",
               "YlGnBu_r",
               "YlOrRd",
               "YlOrRd_r"]

# ---- Set up App ----
ext_css = ["https://use.fontawesome.com/releases/v5.8.1/css/all.css"]
app = dash.Dash(__name__,
                external_stylesheets=ext_css,
                url_base_pathname="/mosaica/",
                suppress_callback_exceptions=True)
server = app.server

#
# Layout
# ------------------------------------------------------------------------------

# --- header ---
header = html.Div(className="head", children=[
    html.H1(children=[
            html.Span(className="fas fa-atom"), " ",
            __title__
            ]),
    html.H2(__subtitle__)
])

# --- Footer ---
footer = html.Div(className="foot", children=[
    html.Div(className="container", children=[
        html.Div(className="about", children=[
            html.H5("About:"),
            html.P([html.Span(className="fas fa-user"), " ",
                    html.A("Germain Salvato Vallverdu", 
                           href="https://gsalvatovallverdu.gitlab.io/")]),
            html.P([
                html.A([html.Span(className="fab fa-github"), " @gvallverdu"],
                       href="https://github.com/gVallverdu"),
                " / ",
                html.A([html.Span(className="fab fa-twitter"), " @gvallverdu"],
                       href="https://twitter.com/gvallverdu")
            ]),
        ]),
        html.Div(className="uppa-logo", children=[
            html.A(href="https://www.univ-pau.fr", children=[
                html.Img(
                    src=app.get_asset_url("img/LogoUPPAblanc.png"),
                )
            ])]
        )
    ])
])

# --- Body: main part of the app ---
body = html.Div(className="container", children=[

    # --- store components for the data
    dcc.Store(id="data-storage", storage_type="memory"),

    # -- two sides div
    html.Div(children=[
        # --- dash bio Molecule 3D Viewer
        html.Div(id="dash-bio-container", children=[
            html.H4("Structure"),

            # --- controls
            html.Div(className="control-panel", children=[
                # --- upload
                html.Div(className="control-label",
                         children="Upload xyz file"),
                dcc.Upload(
                    id='file-upload',
                    children=html.Div(
                        className="upload-area control",
                        children="Drag and Drop or click to select file"
                    ),
                ),

                # --- select data to plot
                html.Div(className="control-label", children="Select data"),
                dcc.Dropdown(
                    className="control",
                    id='dropdown-data',
                    placeholder="Select data"
                ),

                # --- select colorscale
                html.Div(className="control-label",
                         children="Select colormap"),
                dcc.Dropdown(
                    className="control",
                    id='dropdown-colormap',
                    options=[{"label": cm, "value": cm}
                             for cm in plt.cm.cmap_d],
                    value="cividis"
                ),
                html.P("Click on atoms to highlight the corresponding lines"
                       " in the table on the right."),
            ]),
            html.Div(id="dash-bio-viewer"),
            dcc.Graph(id='colorbar', config=dict(displayModeBar=False))
        ]),

        # --- Data table
        html.Div(id="data-table-container", children=[
            html.H4("Data Table"),
            html.Div(className="control-panel", children=[
                html.Div(className="column-selector-label",
                         children="Select the columns of the table:"),
                dcc.Dropdown(
                    id="data-column-selector",
                    multi=True,
                ),
            ]),
            html.Div(children=[
                dash_table.DataTable(
                    id="data-table",
                    editable=False,
                    style_cell={'maxWidth': '60px',
                                'width': '60px',
                                'minWidth': '60px',
                                'whiteSpace': 'normal'},
                    style_header={'backgroundColor': 'rgba(60, 93, 130, .25)',
                                  'fontWeight': 'bold',
                                  "border": "1px solid gray",
                                  "fontFamily": "sans-serif"},
                    style_data_conditional=[{'if': {'row_index': 'odd'},
                                             'backgroundColor': 'rgba(60, 93, 130, .05)'}],
                    style_data={'border': '1px solid gray'},
                    style_table={"overflowX": "scroll",
                                 "maxHeight": "800px",
                                 "overflowY": "scroll"},
                    # need higher version of dash_table, incompatible with dash-bio
                    # fixed_rows={'headers': True, 'data': 0}
                )
            ]),
        ])
    ]),
    html.Hr(className="clearfix"),
    html.Div(className="documentation", children=[
        dcc.Markdown(__doc__)
    ])
])

app.layout = html.Div([header, body, footer])

#
# callbacks
# ------------------------------------------------------------------------------


@app.callback(
    [Output("data-storage", "data"),
     Output("dash-bio-viewer", "children"),
     Output("dropdown-data", "options"),
     Output("data-column-selector", "options"),
     Output("data-column-selector", "value")],
    [Input("file-upload", "contents")]
)
def upload_data(content):
    """
    Uploads the data from an xyz file and store them in the store component.
    Then set up the dropdowns, the table and the molecule viewer.
    """

    # read file
    if content:
        content_type, content_str = content.split(",")
        decoded = base64.b64decode(content_str).decode("utf-8")
        fdata = io.StringIO(decoded)
        species, coords = utils.read_molecule(fdata)

    else:
        filename = app.get_asset_url("data/C28-D2.xyz")
        with open(filename, "r") as f:
            species, coords = utils.read_molecule(f)

    # comute data
    df, distances = utils.compute_data(species, coords)
    model_data = utils.get_molecular_data(species, coords)

    # all data for the store component
    all_data = df.to_dict("records")

    # Set the molecule 3D Viewer component
    dbviewer = dash_bio.Molecule3dViewer(
        id='molecule-viewer',
        backgroundColor="#FFFFFF",
        # backgroundOpacity='0',
        modelData=model_data,
        atomLabelsShown=True,
        selectionType='atom'
    )

    # dropdown for table columns
    tab_options = [{"label": name, "value": name} for name in df]
    values = ["atom index", "species", "angular defect", "haddon", "neighbors"]

    # options for dropdown for mapped values
    options = [{"label": name, "value": name} for name in df
               if name not in ["atom index", "species"]]

    return all_data, dbviewer, options, tab_options, values


@app.callback(
    [Output("data-table", "data"),
     Output("data-table", "columns")],
    [Input("data-storage", "modified_timestamp"),
     Input("data-column-selector", "value")],
    [State("data-storage", "data")]
)
def set_table_columns(ts, values, data):
    """
    Select columns for the table
    """

    df = pd.DataFrame(data)

    if values is None:
        return [], []
    else:
        tab_df = df[values]
        columns = [{"name": i, "id": i} for i in tab_df]
        data = tab_df.to_dict("records")

        return data, columns


@app.callback(
    Output("data-table", "style_data_conditional"),
    [Input("molecule-viewer", "selectedAtomIds")]
)
def highlight_selected_atoms(atom_ids):
    """
    Highlights the columns corresponding to the selected atoms.
    """

    style_data_conditional = [{'if': {'row_index': 'odd'},
                               'backgroundColor': 'rgba(60, 93, 130, .05)'}]
    if atom_ids:
        for iat in atom_ids:
            style_data_conditional.append({
                "if": {"row_index": iat},
                "backgroundColor": 'rgba(60, 93, 130, .75)',
                "color": "white",
            })

    return style_data_conditional


@app.callback(
    Output('molecule-viewer', 'styles'),
    [Input('dropdown-data', 'value'),
     Input('dropdown-colormap', "value")],
    [State("data-storage", "data")]
)
def map_data_on_atoms(selected_data, cm_name, data):
    """
    Map the selected data on the structure using a colorscale to draw the atoms.
    """

    df = pd.DataFrame(data)

    if selected_data:
        normalize = mpl.colors.Normalize(df[selected_data].min(),
                                         df[selected_data].max())
        cm = plt.cm.get_cmap(cm_name)
        norm_cm = cm(X=normalize(df[selected_data].values), alpha=1)
        colors = [mpl.colors.rgb2hex(color) for color in norm_cm]
        styles_data = {
            str(iat): {
                "color": colors[iat],
                "visualization_type": "stick"
            }
            for iat in range(len(df))
        }

    else:
        styles_data = {
            str(iat): {
                "color": utils.get_atom_color(df.species[iat]),
                "visualization_type": "stick"
            }
            for iat in range(len(df))
        }

    return styles_data


@app.callback(
    Output("colorbar", "figure"),
    [Input('dropdown-data', 'value'),
     Input('dropdown-colormap', 'value')],
    [State("data-storage", "data")]
)
def plot_colorbar(selected_data, cm_name, data):
    """
    Display a colorbar according to the selected data mapped on to the structure.
    """

    if selected_data:
        # get data and boundaries
        values = pd.DataFrame(data)[selected_data].values
        minval, maxval = np.nanmin(values), np.nanmax(values)

        # set up fake data and compute corresponding colors
        npts = 100
        values = np.linspace(minval, maxval, npts)
        normalize = mpl.colors.Normalize(minval, maxval)

        cm = plt.cm.get_cmap(cm_name)
        cm_RGBA = cm(X=normalize(values), alpha=1) * 255
        cm_rgb = ["rgb(%d, %d, %d)" % (int(r), int(g), int(b))
                  for r, g, b, a in cm_RGBA]
        colors = [[x, c] for x, c in zip(np.linspace(0, 1, npts), cm_rgb)]

        trace = [
            go.Contour(
                z=[values, values],
                x0=values.min(),
                dx=(values.max() - values.min()) / npts,
                colorscale=colors,
                autocontour=False,
                showscale=False,
                contours=go.contour.Contours(coloring="heatmap"),
                line=go.contour.Line(width=0),
                hoverinfo="skip",
            ),
        ]
        figure = go.Figure(
            data=trace,
            layout=go.Layout(
                width=550, height=100,
                xaxis=dict(showgrid=False, title=selected_data),
                yaxis=dict(ticks="", showticklabels=False),
                margin=dict(l=0, t=0, b=40, r=0, pad=0)
            )
        )
    else:
        figure = go.Figure(
            data=[],
            layout=go.Layout(
                width=550, height=100,
                xaxis=dict(ticks="", showticklabels=False, showgrid=False,
                           title=selected_data, zeroline=False),
                yaxis=dict(ticks="", showticklabels=False, showgrid=False,
                           title=selected_data, zeroline=False),
                margin=dict(l=5, t=0, b=40, r=5, pad=0)
            )
        )

    return figure


if __name__ == '__main__':
    app.run_server(debug=True)
