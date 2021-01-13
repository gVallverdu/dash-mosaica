#!/usr/bin/env python3
# -*- coding=utf-8 -*-


import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
import dash_bio


model_data = {
    'atoms': [
        {'name': 'H1', 'serial': 0, 'element': 'H',
         'positions': [ 0.      ,  0.756697, -0.52061 ]},
        {'name': 'O2', 'serial': 1, 'element': 'O',
         'positions': [0., 0., 0.]},
        {'name': 'H3', 'serial': 2, 'element': 'H',
         'positions': [ 0.      , -0.756697, -0.52061 ]}
    ],
    'bonds': [
        {'atom1_index': 0, 'atom2_index': 1},
        {'atom1_index': 1, 'atom2_index': 2}
    ]
}

# ---- Set up App ----

app = dash.Dash(__name__)

#
# Layout
# ------------------------------------------------------------------------------

# --- Body: main part of the app ---
app.layout = html.Div(className="container", children=[

    # --- store components for the data
    dcc.Store(id="data-storage", storage_type="memory"),

    # -- two sides div
    html.Div(children=[
        # --- dash bio Molecule 3D Viewer
        html.Div(id="dash-bio-container", children=[

            html.Div(id="version"),

            # --- controls
            html.Button("Update", id="button"),

            html.Div(
                dash_bio.Molecule3dViewer(
                    id='molecule-viewer',
                    modelData={"atoms": [], "bonds": []},
                )
                # id="dash-bio-viewer"
            ),
        ]),
    ]),
])


#
# callbacks
# ------------------------------------------------------------------------------

@app.callback(
    Output("version", "children"),
    [Input("button", "n_clicks")]
)
def version(n_clicks):
    line = "Modules' version\n"
    if n_clicks:
        for module in [dash, dash_bio, html, dcc]:
            line += "%s: %s\n" % (module.__version__, module.__name__)
    return line
        

@app.callback(
    Output("molecule-viewer", "modelData"),
    [Input("button", "n_clicks")]
)
def upload_data(n_clicks):
    """
    Uploads the data from an xyz file and store them in the store component.
    """

    # read file
    if n_clicks:
        return model_data
    else:
        return {"atoms": [], "bonds": []}


if __name__ == '__main__':
    app.run_server(debug=True)
