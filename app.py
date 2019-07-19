#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
## Documentation

This application aims to visualize local geometric informations about a
molecular structure. In particular, the application computes geometric
quantities which provide an insight of the discrete curvature of a molecular
structure. In addition, from upload or by editing the table, you can visualize 
any atomic properties.

### Global overview

The tools on the right control the visualization. The options on the left
allow to select the column displayed in the table.

#### On the left

* The dashed box allows to upload the xyz file.
* The *"Select data"* menu, allows to select the data you want to visualize on 
the structure.
* The *"Select colormap"* menu, changes the colormap. The `_r` label 
corresponds to colormap in inverse order.
* The *"colormap boundaries"* inputs change the min and max values used to 
compute the colors associated to the data.

#### On the right

The right panel displays a table of the data. Select the columns you want to 
show using the check boxes. The value in the table can be modified and the 
visualization is updated each time you modify a value.

If you want to add manualy custom data, you can add the `custom` column to the
table and fill it with your values. You can copy and pasta data from a 
spreadsheet or a text file.

The whole data can be downloaded in csv format from the button at the top.

**Warning:** If you edit the data in the table, you have first to refresh the
application before uploading a new molecule.

### Geometrical data

The definitions of the geometrical data available by default are given below.
Some of them are available only if there is a minimum number of bonds:

* **Angular defect (degrees):** The angular defect is a measure of the discrete curvature
on a given atom. It is computed as360° minus the sum of the angles between bonds with atoms
bonded to the considered atom.
* **haddon (degrees):** This is the pyramidalization angle as defined by 
[R.C. Haddon, _C60: Sphere or Polyhedron?_, J. Am. Chem. Soc. **1997**](https://pubs.acs.org/doi/10.1021/ja9637659)
* **improper angle (degrees):** This is the improper dihedral angle
* **dist. from. ave. plane (angstrom):** This is the distance between the
considered atom and the average plane defined by atoms bonded to it.
* **neighbors (number):** This is the number of neighbors of the atom
* **ave. neighb. dist. (angstrom):** This is the average distance between the 
considered atom and its neighbors.

### File and data upload

The application accepts standard xyz files.
Such a file is suposed to display the number of atoms on the first line,
followed by a title line and followed by the structure in cartesian
coordinates. Each line contains the element as first column and the
cartesian coordinates as 2d, 3th and 4th columns, for example:

    3
    H2O molecule
    O   -0.111056  0.033897  0.043165
    H    0.966057  0.959148 -1.089095
    H    0.796629 -1.497157  0.403985


If additional data are provided on each line, they are returned as
atomic properties and can be visualized. The names of theses additional atomic
properties are `propX`, `X` being a number. For example, the file below will provide
an atomic properties with name `prop0`:

    3
    H2O molecule
    O   -0.111056  0.033897  0.043165   -1.8
    H    0.966057  0.959148 -1.089095    0.9
    H    0.796629 -1.497157  0.403985    0.9

The number of provided data must be the same on each line. If this is not the
case, only the structure is read.

"""

import io
import base64
import urllib

import dash
import dash_table
from dash_table.Format import Format, Scheme
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
    html.H1(children=[html.Span(className="fas fa-atom"), " ", __title__]),
    # html.H2(__subtitle__)
    html.A(
        id="github-link",
        href="https://github.com/gVallverdu/dash-mosaica",
        children=[
            "View on GitHub",
        ]
    ),
    html.Span(id="github-icon", className="fab fa-github fa-2x"),
])

# --- Footer ---
footer = html.Div(className="foot", children=[
    html.Div(className="container", children=[
        html.Div(className="about", children=[
            html.H5("About:"),
            html.P([
                html.A("Germain Salvato Vallverdu",
                       href="https://gsalvatovallverdu.gitlab.io/")]),
            html.P(
                html.A(href="https://www.univ-pau.fr", children=[
                    "University of Pau & Pays Adour"
                ])
            )
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
                html.Div(id="upload-label", children="Upload xyz file"),
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

                # --- select colormap
                html.Div(className="control-label",
                         children="Select colormap"),
                dcc.Dropdown(
                    className="control",
                    id='dropdown-colormap',
                    options=[{"label": cm, "value": cm}
                             for cm in plt.cm.cmap_d],
                    value="cividis"
                ),

                # --- colormap boundaries
                html.Div(className="control-label",
                         children="Colormap boundaries"),
                html.Div(className="control", children=[
                    dcc.Input(id="cm-min-value", type="number", debounce=True,
                              placeholder=0),
                    dcc.Input(id="cm-max-value", type="number", debounce=True,
                              placeholder=0),
                ]),

                html.P("Click on atoms to highlight the corresponding lines"
                       " in the table on the right."),
            ]),
            dcc.Graph(id='colorbar', config=dict(displayModeBar=False)),
            html.Div(id="dash-bio-viewer"),
        ]),

        # --- Data table
        html.Div(id="data-table-container", children=[
            html.A(
                html.Button("Download data", id="download-button"),
                id="download",
                download="rawdata.csv",  href="",
                target="_blank",
            ),
            html.H4("Data Table"),
            html.Div(className="control-panel", children=[
                html.Div(className="column-selector-label",
                         children="Select the columns of the table:"),
                dcc.Checklist(
                    id="data-column-selector",
                    value=[],
                    inputClassName="checklist-item",
                    labelClassName="checklist-label",
                ),
            ]),
            html.Div(children=[
                dash_table.DataTable(
                    id="data-table",
                    editable=True,
                    row_selectable="multi",
                    style_cell={'minWidth': '60px', 'whiteSpace': 'normal'},
                    style_header={
                        'backgroundColor': 'white',
                        "padding": "5px",
                        'fontWeight': 'bold',
                        "textAlign": "center",
                        "borderBottom": "2px solid rgb(60, 93, 130)",
                        "borderTop": "2px solid rgb(60, 93, 130)",
                        "fontFamily": "sans-serif"},
                    style_data_conditional=[{
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgba(60, 93, 130, .05)'
                    }],
                    style_table={"overflowX": "scroll",
                                 "maxHeight": "700px",
                                 "overflowY": "scroll"},
                    fixed_rows={'headers': True, 'data': 0}
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
    [Input("file-upload", "contents"),
     Input('data-table', 'data_timestamp')],
    [State("data-storage", "data"),
     State("data-table", "data"),
     State("data-column-selector", "value"),
     State("dash-bio-viewer", "children")
     ]
)
def upload_data(content, table_ts, stored_data, table_data, selected_columns,
                dbviewer):
    """
    Uploads the data from an xyz file and store them in the store component.
    Then set up the dropdowns, the table and the molecule viewer.
    """

    if table_ts is not None:
        # update stored data from current data in the table
        df = pd.DataFrame(stored_data)
        try:
            table_df = pd.DataFrame(table_data)
            table_df = table_df.astype({col: np.float for col in table_df
                                        if col != "species"})
            df.update(table_df)
        except ValueError:
            print("No update of data")

        all_data = df.to_dict("records")

    else:
        # Initial set up, read data from upload

        # read file
        if content:
            content_type, content_str = content.split(",")
            decoded = base64.b64decode(content_str).decode("utf-8")
            fdata = io.StringIO(decoded)
            species, coords, atomic_prop = utils.read_xyz(fdata)

        else:
            # filename = app.get_asset_url("data/C28-D2.xyz")
            filename = "assets/data/C28-D2.xyz"
            with open(filename, "r") as f:
                species, coords, atomic_prop = utils.read_xyz(f)

        # comute data
        df, distances = utils.compute_data(species, coords)
        if atomic_prop:
            df_prop = pd.DataFrame(atomic_prop, index=df.index)
            df = pd.merge(df, df_prop, left_index=True, right_index=True)

        if "custom" not in df:
            df["custom"] = 0.0

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

        # options for the checklist in order to select the columns of the table
        selected_columns = ["atom index", "species", "angular defect",
                            "haddon", "neighbors"]

    # options to select data mapped on atoms
    options = [{"label": name, "value": name} for name in df
               if name not in ["atom index", "species"]]

    # checklist options to select table columns
    tab_options = [{"label": name, "value": name} for name in df]

    return all_data, dbviewer, options, tab_options, selected_columns


@app.callback(
    [Output("data-table", "data"),
     Output("data-table", "columns")],
    [Input("data-storage", "modified_timestamp"),
     Input("data-column-selector", "value")],
    [State("data-storage", "data")]
)
def select_table_columns(ts, values, data):
    """
    Select columns displayed in the table. A custom column is available and 
    filled with zero by default.
    """

    # get data from the Store component
    df = pd.DataFrame(data)

    if values is None:
        # initial set up
        return [], []
    else:
        # fill the table with the selected columns
        tab_df = df[values]
        data = tab_df.to_dict("records")

        # add format
        columns = list()
        for column in tab_df:
            if column in {"atom index", "species", "neighbors", "custom"}:
                columns.append({"name": column, "id": column})
            elif column[:4] == "prop":
                columns.append({"name": column, "id": column})
            else:
                columns.append({
                    "name": column, "id": column, "type": "numeric",
                    "format": Format(
                        precision=4,
                        scheme=Scheme.fixed,
                    )
                })

        return data, columns


@app.callback(
    [Output("data-table", "style_data_conditional"),
     Output("data-table", "selected_rows")],
    [Input("molecule-viewer", "selectedAtomIds")]
)
def select_rows_from_atoms(atom_ids):
    """
    Highlights the rows corresponding to the selected atoms in the molecule
    viewer.
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

    return style_data_conditional, list(atom_ids)


# @app.callback(
#     Output("molecule-viewer", "selectedAtomIds"),
#     [Input("data-table", "selected_rows")]
# )
# def select_atom_from_table(rows):
#     """
#     Select atom from rows selected in the table
#     """
#     if rows:
#         print(rows)
#         return list(rows)
#     else:
#         return []


@app.callback(
    Output('molecule-viewer', 'styles'),
    [Input('dropdown-data', 'value'),
     Input('dropdown-colormap', "value"),
     Input("data-storage", "modified_timestamp"),
     Input("cm-min-value", "value"),
     Input("cm-max-value", "value")],
    [State("data-storage", "data")]
)
def map_data_on_atoms(selected_data, cm_name, ts, cm_min, cm_max, data):
    """
    Map the selected data on the structure using a colormap.
    """

    df = pd.DataFrame(data)

    if selected_data:
        values = df[selected_data].values
        minval, maxval = np.nanmin(values), np.nanmax(values)

        # get cm boundaries values from inputs if they exist
        if cm_min:
            minval = cm_min
        if cm_max:
            maxval = cm_max

        normalize = mpl.colors.Normalize(minval, maxval)

        cm = plt.cm.get_cmap(cm_name)
        norm_cm = cm(X=normalize(values), alpha=1)
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
     Input('dropdown-colormap', 'value'),
     Input("data-storage", "modified_timestamp"),
     Input("cm-min-value", "value"),
     Input("cm-max-value", "value")],
    [State("data-storage", "data")]
)
def plot_colorbar(selected_data, cm_name, data_ts, cm_min, cm_max, data):
    """
    Display a colorbar according to the selected data mapped on to the structure.
    """

    if selected_data:
        # get data and boundaries
        values = pd.DataFrame(data)[selected_data].values
        minval, maxval = np.nanmin(values), np.nanmax(values)

        # get cm boundaries values from inputs if they exist
        if cm_min:
            minval = cm_min

        if cm_max:
            maxval = cm_max

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
                dx=(values.max() - values.min()) / (npts - 1),
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
                width=580, height=100,
                xaxis=dict(showgrid=False, title=selected_data),
                yaxis=dict(ticks="", showticklabels=False),
                margin=dict(l=40, t=0, b=40, r=20, pad=0)
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


@app.callback(Output('download', 'href'),
              [Input('data-storage', 'modified_timestamp')],
              [State("data-storage", "data")])
def update_download_button(ts, data):
    """
    Return a link with the data in csv format.
    """

    if ts is not None:
        df = pd.DataFrame(data)

        # put atom index, species, and coordinates as first columns
        first_columns = ["atom index", "species", "x", "y", "z"]
        columns = df.columns.drop(first_columns)
        df = df[first_columns + list(columns)]

        csv_string = df.to_csv(index=False, encoding="utf-8")
        csv_string = "data:text/csv;charset=utf-8," + \
            urllib.parse.quote(csv_string)
        return csv_string

    else:
        return "#"


if __name__ == '__main__':
    app.run_server(debug=True)
