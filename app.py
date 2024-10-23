#Module import
from dash import Dash, html, dcc, Output, Input, State, ctx, dash_table
import dash_daq as daq
import plotly.graph_objects as go
import skywalker
import pandas as pd
from datetime import date
from astropy.time import Time

#=--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--=

# Initialize the app
#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
#app = Dash(__name__, external_stylesheets=external_stylesheets)
app = Dash(__name__)
server = app.server

#Creating night instance
n = skywalker.night('aao')

#Creating the title and description card
def title_card():
    """
    :return: A Div containing dashboard title & descriptions.
    """
    return html.H1("SkyWalker", style={'font-family':'Trebuchet MS','margin-top':20,'margin-left':20, 'color':'#ff9f06', 'font-weight':'bold'})

#Creating the title and description card
def description_card():
    """
    :return: A Div containing dashboard title & descriptions.
    """
    return html.Div(
        id="description-card",
        children=[
            html.Div(
                id="intro",
                children=[
                    "SkyWalker was written to enable visualisation of target observability for astronomical observations, ",
                    "as well as generate optimised observing queues taking into account limiting elevation angles, lunar angular distances, and pointing restrictions. ",
                    "Upload a .csv file with columns as 'name', 'ra', 'dec' to get started. ",
                    "Columns of 'priority' and 'obs_time' are also required for queue generation."
                    ],
                style={'font-family':'Trebuchet MS','margin-left':20, 'margin-right':20,'margin-bottom':20,'margin-top':20, 'color':'#DDDDDD'}
                ),
            html.Div(
                id="control_instructions_header",
                children=[
                    "The control panel functions as follows:"
                    ],
                style={'font-family':'Trebuchet MS','margin-left':20, 'margin-right':20,'margin-bottom':0, 'margin-top':40, 'color':'#DDDDDD'}
                ),
            html.Div(
                id="control_instructions",
                children=[
                    "-Date: select the UT date when the local night begins.", html.Br(),
                    "-Observatory: Select the observatory from the provided list (the introduction of custom coordinates to be implemented).", html.Br(),
                    "-Limiting twilight: Select the twlight threshold from which to begin and conclude the queue generation.", html.Br(),
                    "-Minimum/maximum observing angle: Select the minimum and maximum elevation angle thresholds for queue generation.", html.Br(),
                    "-Minimum lunar distance: Select the minimum allowed angular distance from the moon (in degrees) for queue generation.", html.Br(),
                    "-Unallowed pointings: Select nautical directions in which targets should not be observed for queue generation.", html.Br(),
                    ],
                style={'font-family':'Trebuchet MS','margin-left':40, 'margin-right':20,'margin-bottom':20, 'margin-top':0, 'color':'#DDDDDD'}
                )
            ]
        )

#Toggle lunar distances
#daq.BooleanSwitch(
#    id = 'lunar_distance_toggle',
#    label = 'Lunar distances',
#    labelPosition = 'top',
#    on = False,
#    color = '#ae93ff',
#    style={"color": "#dddddd", 'margin-bottom':5}),

#Creating the upload card
def upload_card():
    """
    :return: A Div containing dashboard controls.
    """
    return html.Div(
        children=[
            #Text input
            html.P('Target entry',
                style={
                    'textAlign':'center',
                    'font-size':14,
                    'color':'#dddddd',
                    'margin-top':10
                    }
                ),
            html.Div(
                dcc.Textarea(
                    id='add_target_text',
                    value='',
                    placeholder = 'ra dec name priority time',
                    style={'width': '100%', 'height': 100},
                    ),
                style = {'margin-left':'5%', 'margin-right':'5%'}
                ),
            html.Center(html.Button('Add', id='add_target_button', n_clicks=0)),

            #File upload
            html.P('File upload',
                style={
                    'textAlign':'center',
                    'font-size':14,
                    'color':'#dddddd',
                    'margin-top':30
                    }
                ),
            dcc.Upload(
                id="upload_data",
                children=html.Div(['Drag and Drop or ', html.A('Select File')]),
                style={
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "margin": "10px",
                    "color":"#DDDDDD"
                    },
                multiple=False
                ),

            #Clear targets
            html.Center(
                html.Button('Clear targets', id='clear_target_button', n_clicks=0),
                style = {'margin-top':40})
            ]
            )

#Creating the controls card
def controls_card():
    """
    :return: A Div containing dashboard controls.
    """
    return html.Div(
        children=[
            #Date selection
            html.P('Date',
                style={
                    'textAlign':'center',
                    'font-size':14,
                    'color':'#dddddd',
                    'margin-top':10
                    }
                ),
            html.Center([
                dcc.DatePickerSingle(
                    id='date_picker',
                    min_date_allowed=date(2000, 1, 1),
                    max_date_allowed=date(2080, 12, 31),
                    initial_visible_month=date.today(),
                    date=date.today(),
                    display_format='DD MMMM Y',
                    first_day_of_week = 1
                    )],
                style = {'margin-bottom':5}
                ),

            #Observatory selection
            html.P('Observatory',
                style={
                    'textAlign':'center',
                    'font-size':14,
                    'color':'#dddddd'
                    }
                ),
            html.Center(
                dcc.Dropdown(
                    options=list(observatories.index),
                    value = 'Australian Astronomical Observatory',
                    id='observatory_dropdown',
                    clearable=False),
                    style = {'margin-left': '5%', 'margin-right': '5%', 'margin-bottom':5}
                ),

            #Limiting twilight selection
            html.P('Limiting twilight',
                style={
                    'textAlign':'center',
                    'font-size':14,
                    'color':'#dddddd'
                    }
                ),
            html.Center(
                dcc.Dropdown(
                    options=['Astronomical','Nautical','Civil'],
                    value = 'Astronomical',
                    id='twilight_dropdown',
                    clearable=False),
                    style = {'margin-left': '5%', 'margin-right': '5%', 'margin-bottom':5}
                ),

            #Threshold angle selection
            html.Center(
                daq.NumericInput(
                    id='minimum_angle',
                    value=0,
                    min = 0,
                    max = 89,
                    label = 'Minimum observation angle',
                    style={"color": "#dddddd"},
                    className='NumericInput'
                    ),
                style = {'margin-bottom':5}
                ),
            html.Center(
                daq.NumericInput(
                    id='maximum_angle',
                    value=90,
                    min = 1,
                    max = 90,
                    label = 'Maximum observation angle',
                    style={"color": "#dddddd"},
                    className='NumericInput'
                    ),
                style = {'margin-bottom':5}
                ),

            #Setting lunar distance threshold
            html.Center(
                daq.NumericInput(
                    id='minimum_lunar_distance',
                    value=20,
                    min = 0,
                    max = 359,
                    label = 'Minimum lunar distance',
                    style={"color": "#dddddd"},
                    className='NumericInput'
                    ),
                style = {'margin-bottom':5}
                ),
            html.P(
                'Unallowed pointings',
                style={
                    'textAlign':'center',
                    'font-size':14,
                    'color':'#dddddd'
                    }
                ),

            #Selecting unallowed pointings
            html.Center([dcc.Checklist(
                            ['N','NE','E','SE','S','SW','W','NW'],
                            value = [],
                            id = 'pointing_checklist',
                            inline = True,
                            style = {'align':'center', 'margin-bottom':20, 'color':'#dddddd'})]),
            html.Center(html.Button('Generate queue', id='gen_queue', n_clicks=0))
            ]
        )

#Compiling the card components to the overall application layout
def layout():
    """
    :return: List of layout components.
    """
    return [
        #Left column
        html.Div(
            id="left-column",
            className="three columns",
            children = [
                title_card(),
                dcc.Tabs([
                    dcc.Tab(
                        label = 'Upload',
                        className='custom-tab-light',
                        selected_className='custom-tab-light-selected',
                        children = [upload_card()]),
                    dcc.Tab(
                        label = 'Controls',
                        className='custom-tab-light',
                        selected_className='custom-tab-light-selected',
                        children = [controls_card()])
                    ])
                ],
            style={'background-color': '#2b2b2b', 'font-family':'Trebuchet MS', "height": "100vh"}
            ),
        #Right column
        html.Div(
            id="right-column",
            className="eight columns",
            children=[
                dcc.Tabs([
                    #Observability plot tab
                    dcc.Tab(
                        label = 'Observability plot',
                        className='custom-tab',
                        selected_className='custom-tab-selected',
                        children = [dcc.Graph(id = 'graph_paths', figure = fig_paths, style={'height':600})]),

                    #Targets table tab
                    dcc.Tab(
                        label = 'Targets',
                        className='custom-tab',
                        selected_className='custom-tab-selected',
                        children = [
                            html.Div(
                                dash_table.DataTable(
                                    id='data_table',
                                    data=n.data.to_dict('records'),
                                    columns=[{"name": 'Name', "id": 'name'}, {"name": 'RA', "id": 'ra'}, {"name": 'DEC', "id": 'dec'}],
                                    export_columns='visible',
                                    export_format='csv',
                                    style_header={
                                        'backgroundColor': '#1e1e1e',
                                        'color': '#ddd',
                                        'textAlign': 'center'
                                        },
                                    style_data={
                                        'backgroundColor': '#1e1e1e',
                                        'color': '#ddd',
                                        'textAlign': 'left'
                                        }
                                ), 
                            style = {'margin-top':20}
                            )]
                        ),

                    #Queue plot tab
                    dcc.Tab(
                        label = 'Queue plot',
                        className='custom-tab',
                        selected_className='custom-tab-selected',
                        children = [
                            dcc.Graph(id = 'graph_queue', figure = fig_queue, style={'height':600})
                        ]),

                    #Queue table tab
                    dcc.Tab(
                        label = 'Queue',
                        className='custom-tab',
                        selected_className='custom-tab-selected',
                        children = [
                            html.Div(
                                dash_table.DataTable(
                                    id='queue_table',
                                    data=n.queue.to_dict('records'),
                                    columns=[{"name": 'Name', "id": 'name'}, {"name": 'Priority', "id": 'priority'}, {"name": 'RA', "id": 'ra'}, {"name": 'DEC', "id": 'dec'}, {"name": 'Observation time', "id": 'obs_times'}, {"name": 'UT start time', "id": 'start_ut'}],
                                    export_columns='visible',
                                    export_format='csv',
                                    style_header={
                                        'backgroundColor': '#1e1e1e',
                                        'color': '#ddd',
                                        'textAlign': 'center'
                                        },
                                    style_data={
                                        'backgroundColor': '#1e1e1e',
                                        'color': '#ddd',
                                        'textAlign': 'left'
                                        }
                                    ),
                                style = {'margin-top':20}
                                )
                        ]),

                    #Instructions tab
                    dcc.Tab(
                        label = 'Instructions',
                        className='custom-tab',
                        selected_className='custom-tab-selected',
                        children = [description_card()])
                    ], style={'margin-top':100})
                ],
            style={'background-color': '#1e1e1e', 'font-family':'Trebuchet MS', "height": "100vh"}
            ),
        dcc.Store(id = 'night-object')
        ]

fig_paths = n.plot_paths(dashapp = True, dark = True)
fig_queue = n.plot_paths(dashapp = True, dark = True)

#Importing data for observatories dropdown selection
observatories = pd.read_csv('./observatories.csv', delimiter = '\t', index_col = 'Name')

#Defining app_layout
app.layout = layout()

#=--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--=

#User input observability plot and data table tabs
@app.callback(
    Output("graph_paths", "figure"),
    Output("data_table", "data"),
    Output("upload_data", "filename"),
    Output("upload_data", "contents"),
    Output("add_target_text", "value"),
    Output("night-object", "data"),

    Input("night-object", "data"),
    Input("date_picker", "date"),
    Input("observatory_dropdown", "value"),
    Input("upload_data", "contents"),
    Input("minimum_angle", "value"),
    Input("maximum_angle", "value"),
    Input("add_target_button", 'n_clicks'),
    State('add_target_text', 'value'),
    Input("minimum_lunar_distance", "value"),
    Input("twilight_dropdown", "value"),
    Input("pointing_checklist", 'value'),
    Input("clear_target_button", 'n_clicks')
    )
def update_paths_plot(store_input, date, observatory, upload_data, min_angle, max_angle, add_target, add_target_text, minimum_lunar_distance, limiting_twilight, restricted_directions, clear_target):  

    n = skywalker.night(
        obs_id = observatories['ID'][observatory], 
        date = date, 
        min_angle = min_angle, 
        max_angle = max_angle, 
        minimum_lunar_distance = minimum_lunar_distance, 
        limiting_twilight = limiting_twilight.lower(), 
        restricted_directions = restricted_directions)

    if store_input != None:
        if len(store_input['names']) != 0:
            n.load_targets(objects = store_input['target_coords'], names = store_input['names'], obs_times = store_input['obs_times'], priorities = store_input['priorities'])

    """
    #Updating date
    if date != n.date_ymd:
        n.change_date(date)
    #Updating observatory
    if observatory != n.observatory.name:
        n.change_observatory(observatories['ID'][observatory])
    #Updating minimum angle
    if min_angle != n.min_angle:
        n.change_limiting_angle('min', min_angle)
    #Updating maximum angle
    if max_angle != n.max_angle:
        n.change_limiting_angle('max', max_angle)
    #Updating minimum lunar angular distance
    if minimum_lunar_distance != n.minimum_lunar_distance:
        n.change_limiting_angle('moon', minimum_lunar_distance)
    #Updating limiting twilight
    if limiting_twilight != n.limiting_twilight:
        n.change_limiting_angle('twilight', limiting_twilight)
    #Updating unallowed pointings
    if restricted_directions != n.restricted_directions:
        n.change_restricted_directions(restricted_directions)
    """

    #Uploading targets
    if upload_data != None:
        n.load_targets_dashapp(upload_data)

    if "clear_target_button" == ctx.triggered_id:
        n.clear_targets()

    if "add_target_button" == ctx.triggered_id:
        input_lines = add_target_text.split(sep = '\n')

        for t in input_lines:
            inputs = t.split(sep = ' ')
    
            try:
                ra = float(inputs[0])
                dec = float(inputs[1])
            
                if len(inputs) > 2:
                    name = inputs[2]
                else:
                    name = None
    
                if len(inputs) > 3:
                    try:
                        priority = int(inputs[3])
                    except:
                        priority = None
                else:
                    priority = None
    
                if len(inputs) > 4:
                    try:
                        obs_time = int(inputs[4])
                    except:
                        obs_time = None
                else:
                    obs_time = None
            
                n.add_target(ra, dec, name = name, obs_time = obs_time, priority = priority)

            except:
                pass

        output_text = ''
    else:
        output_text = add_target_text

    #Updating output table
    output_table = n.data.copy().reset_index()

    #Creating store output dictionary
    store_output = {
        'target_coords':n.target_coords,
        'names':n.names,
        'priorities':n.priorities,
        'obs_times':n.obs_times,
        }

    return n.plot_paths(dashapp = True, dark = True), output_table.to_dict('records'), '', None, output_text, store_output

#User input observability plot and data table tabs
@app.callback(
    Output("graph_queue", "figure"),
    Output("queue_table", "data"),

    Input("night-object", "data"),
    Input("gen_queue", 'n_clicks'),
    Input("date_picker", "date"),
    Input("observatory_dropdown", "value"),
    Input("upload_data", "contents"),
    Input("minimum_angle", "value"),
    Input("maximum_angle", "value"),
    Input("add_target_button", 'n_clicks'),
    Input('add_target_text', 'value'),
    Input("minimum_lunar_distance", "value"),
    Input("twilight_dropdown", "value"),
    Input("pointing_checklist", 'value'),
    Input("clear_target_button", 'n_clicks')
    )
def update_queue_plot(store_input, trigger_queue, date, observatory, upload_data, min_angle, max_angle, add_target, add_target_text, minimum_lunar_distance, limiting_twilight, restricted_directions, clear_target):

    n = skywalker.night(
        obs_id = observatories['ID'][observatory], 
        date = date, 
        min_angle = min_angle, 
        max_angle = max_angle, 
        minimum_lunar_distance = minimum_lunar_distance, 
        limiting_twilight = limiting_twilight.lower(), 
        restricted_directions = restricted_directions)

    if store_input != None:
        if len(store_input['names']) != 0:
            n.load_targets(objects = store_input['target_coords'], names = store_input['names'], obs_times = store_input['obs_times'], priorities = store_input['priorities'])

    if "gen_queue" == ctx.triggered_id:
        n.generate_queue(dashapp = True)

    output_table = n.queue.copy().reset_index()

    return n.plot_queue(dashapp = True, dark = True), output_table.to_dict('records')

# Run the app
if __name__ == '__main__':
    app.run(debug=True)
