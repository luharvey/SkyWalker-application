# █▀ █▄▀ █▄█ █░█░█ ▄▀█ █░░ █▄▀ █▀▀ █▀█
# ▄█ █░█ ░█░ ▀▄▀▄▀ █▀█ █▄▄ █░█ ██▄ █▀▄

#Module import
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import rc
import pandas as pd
from datetime import date
from skyfield.api import load
from astropy.time import Time
import functions
#from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import duckdb
import plotly.graph_objects as go
import base64
import io

#Disabling warning regarding the setting of values in dataframe copies
pd.options.mode.chained_assignment = None

#Initialising time properties
ts = load.timescale()
today = str(date.today()) + 'T00:00:00.0'
t = Time(today, format = 'isot', scale = 'utc')
today_mjd = t.mjd

#Plot styles
linestyles = {0:'-', 1:'-.', 2:':', 3:'--'}
colours = ['C3', '#FC6A0B', '#FCBE0B', '#00E65B', '#4BE1C3', 'C9', '#488DFF', '#B27CFF', 'C6', '#FF51B5']

#Night class
#This is the class that shall be initialised by the user, to which targets can be loaded and from which plots can be generated
class night():
    def __init__(self, obs_id, date = today_mjd, min_angle = 0, max_angle = 90, minimum_lunar_distance = 0):
        #Initialising arrays
        self.target_coords = []
        self.night_peak = []
        self.names = []
        self.priorities = []
        self.obs_times = []
        self.obs_times_angles = []
        self.paths = {}
        self.data = pd.DataFrame(columns = ['name', 'ra', 'dec', 'priority', 'obs_times'])
        self.min_angle = min_angle
        self.max_angle = max_angle
        self.minimum_lunar_distance = minimum_lunar_distance
        self.queue = pd.DataFrame(columns = ['name', 'ra', 'dec', 'priority', 'obs_times', 'obs_time_angles', 'night_peak', 'start_angle', 'start_ut']).set_index('name')
        
        #Setting observatory properties
        self.observatory = functions.retrieve_observatory(obs_id)
        if self.observatory == -1:
            return
        
        #Setting date property
        #Input date refers to the date when the local night starts
        try:
            if type(date) == str:
                if '-' in date and 'T' not in date:
                    date_ext = str(date) + 'T00:00:00.0'
                    t = Time(date_ext, format = 'isot', scale = 'utc')
                    self.date = float(t.mjd)
                else:
                    t = Time(date, format = 'isot', scale = 'utc')
                    self.date = float(t.mjd)
            else:
                self.date = date
        except:
            print('The date parameter must be given in one of the following formats:\n\tmjd\n\tyyyy-mm-dd\n\tyyyy-mm-ddThh:mm:ss.ss')
            return -1
        self.date_ymd = Time(self.date, format = 'mjd', scale = 'utc').isot.split(sep = 'T')[0]

        #Calculating the base earth rotation angle(ERA)
        self.era0 = functions.mjd_to_era(self.date + 1)

        #Calculating the ERA array
        self.era = np.linspace((self.era0-functions.day_rotation/2)+self.observatory.longitude, (self.era0+functions.day_rotation/2)+self.observatory.longitude, 361)

        #Calculating the mjd array
        self.mjd = functions.era_to_mjd(self.era)

        #Calculating the utc array
        self.utc = functions.era_to_utc(self.era)

        #Querying solar coordinates (moving sun)
        solar_ra_ref = [functions.solar_radec_mjd(self.mjd[0])[0], functions.solar_radec_mjd(self.mjd[-1])[0]]
        if solar_ra_ref[-1]-solar_ra_ref[0] < 0:
            solar_ra_ref[-1] += 360
        solar_dec_ref = [functions.solar_radec_mjd(self.mjd[0])[1], functions.solar_radec_mjd(self.mjd[-1])[1]]
        solar_ra = np.interp(self.mjd, [self.mjd[0], self.mjd[-1]], solar_ra_ref)
        solar_dec = np.interp(self.mjd, [self.mjd[0], self.mjd[-1]], solar_dec_ref)
        self.solar_coords = [list(solar_ra), list(solar_dec)]

        #Calculating twilights and solar path
        self.twilights, self.solar_path = functions.twilights(self.solar_coords, self.observatory, self.era)
        
        #Retrieving lunar coordinates and calculating lunar path (moving moon)
        lunar_ra_ref = [functions.lunar_radec_mjd(self.mjd[0])[0], functions.lunar_radec_mjd(self.mjd[-1])[0]]
        if lunar_ra_ref[-1]-lunar_ra_ref[0] < 0:
            lunar_ra_ref[-1] += 360
        lunar_dec_ref = [functions.lunar_radec_mjd(self.mjd[0])[1], functions.lunar_radec_mjd(self.mjd[-1])[1]]
        lunar_ra = np.interp(self.mjd, [self.mjd[0], self.mjd[-1]], lunar_ra_ref)
        lunar_dec = np.interp(self.mjd, [self.mjd[0], self.mjd[-1]], lunar_dec_ref)
        self.lunar_coords = [list(lunar_ra), list(lunar_dec)]
        self.lunar_path = functions.path_era([self.lunar_coords[0], self.lunar_coords[1]], self.observatory, self.era)
    
    #Loading observation targets    
    def load_targets_dashapp(self, contents):
        content_type, content_string = contents.split(',')

        decoded = base64.b64decode(content_string)

        #if '.csv' in filename:
        data = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        
        if len(data.columns) == 1:
            print('Cannot read objects from ' + filename + '. The default delimiter is set as a comma.')
            return -1
        
        self.target_coords = []
        
        for i in range(len(data['name'])):
            ra0 = functions.convert_hms_ra(data['ra'][i])
            dec0 = functions.convert_hms_dec(data['dec'][i])

            if ra0 < self.era[0]:
                night_peak = int(round(ra0 + 360, 0))
            else:
                night_peak = int(round(ra0, 0))

            self.target_coords.append([ra0, dec0])
            self.night_peak.append(night_peak)
            self.names.append(data['name'][i])

            try:
                self.priorities.append(int(data['priority'][i]))
            except:
                pass
            try:
                self.obs_times.append(data['obs_time'][i])
                self.obs_times_angles.append(int(np.ceil(data['obs_time'][i] * functions.degrees_per_second)))
            except:
                pass
        
        #Defining dictionary from which to construct the dataframe
        data_df = {'name':self.names, 'ra':np.array(self.target_coords)[:,0], 'dec':np.array(self.target_coords)[:,1]}
        if len(self.priorities) == len(self.names):
            data_df['priority'] = self.priorities
        if len(self.obs_times) == len(self.names):
            data_df['obs_times'] = self.obs_times
            data_df['obs_time_angles'] = self.obs_times_angles
        if len(self.night_peak) == len(self.names):
            data_df['night_peak'] = self.night_peak
        
        #Constructing dataframe
        self.data_duckdb = pd.DataFrame(data_df)
        self.data = self.data_duckdb.set_index('name')
        
        #Calculating target paths across the sky
        self.paths = {}
        self.lunar_distances = {}
        self.nautical_directions = {}
        for i,n in enumerate(self.names):
            self.paths[n] = functions.path_era([self.data['ra'][i], self.data['dec'][i]], self.observatory, self.era)
            self.lunar_distances[n] = functions.angle_to_moon([self.data['ra'][i], self.data['dec'][i]], self.lunar_coords)
            self.nautical_directions[n] = functions.calculate_pointing([self.data['ra'][i], self.data['dec'][i]], self.observatory, self.era)

    #Loading observation targets
    def load_targets(self, objects, names = None, obs_times = None, priorities = None, delimiter = ','):
        #Parsing of list of coordinates
        if type(objects) == list:
            self.target_coords = objects

            for i in range(len(self.target_coords)):
                if self.target_coords[i][0] < self.era[0]:
                    self.night_peak.append(int(round(self.target_coords[i][0] + 360, 0)))
                else:
                    self.night_peak.append(int(round(self.target_coords[i][0], 0)))
            
            #Names (default to numbered list)
            if type(names) == list:
                if len(names) == len(objects):
                    self.names = names
                else:
                    print('The names and object lists must be equal in length.')
                    return -1
            else:
                self.names = list(np.arange(1, len(objects)+1))

            #Priorities (default to 1)
            if type(priorities) == list:
                if len(priorities) == len(objects):
                    self.priorities = priorities
                else:
                    print('The priorities and object lists must be equal in length.')
                    return -1
            else:
                self.priorities = list(np.ones(len(objects)))

            for i in range(len(self.priorities)):
                self.priorities[i] = int(self.priorities[i])

            #Observation times (default to 600s)
            if type(obs_times) == list:
                if len(obs_times) == len(objects):
                    self.obs_times = obs_times
                    self.obs_times_angles = list(np.array(obs_times) * functions.degrees_per_second)
                else:
                    print('The obs_times and object lists must be equal in length.')
                    return -1
            else:
                self.ob_times = list(np.zeros_like(objects) + 600)
                self.obs_times_angles = list(np.zeros_like(objects) + np.ceil(600 * functions.degrees_per_second))

            for i in range(len(self.obs_times_angles)):
                self.obs_times_angles[i] = int(np.ceil(self.obs_times_angles[i]))
        
        #Parsing of input csv file
        elif type(objects) == str:
            data = pd.read_csv(objects, sep = delimiter)
            
            if len(data.columns) == 1:
                print('Cannot read objects from ' + filename + '. The default delimiter is set as a comma.')
                return -1
            
            self.target_coords = []
            
            for i in range(len(data['name'])):
                ra0 = functions.convert_hms_ra(data['ra'][i])
                dec0 = functions.convert_hms_dec(data['dec'][i])

                if ra0 < self.era[0]:
                    night_peak = int(round(ra0 + 360, 0))
                else:
                    night_peak = int(round(ra0, 0))

                self.target_coords.append([ra0, dec0])
                self.night_peak.append(night_peak)
                self.names.append(data['name'][i])

                try:
                    self.priorities.append(int(data['priority'][i]))
                except:
                    pass
                try:
                    self.obs_times.append(data['obs_time'][i])
                    self.obs_times_angles.append(int(np.ceil(data['obs_time'][i] * functions.degrees_per_second)))
                except:
                    pass
        
        #Defining dictionary from which to construct the dataframe
        data_df = {'name':self.names, 'ra':np.array(self.target_coords)[:,0], 'dec':np.array(self.target_coords)[:,1]}
        if len(self.priorities) == len(self.names):
            data_df['priority'] = self.priorities
        if len(self.obs_times) == len(self.names):
            data_df['obs_times'] = self.obs_times
            data_df['obs_time_angles'] = self.obs_times_angles
        if len(self.night_peak) == len(self.names):
            data_df['night_peak'] = self.night_peak

        #Constructing dataframe
        self.data_duckdb = pd.DataFrame(data_df)
        self.data = self.data_duckdb.set_index('name')
        
        #Calculating target paths across the sky
        self.paths = {}
        self.lunar_distances = {}
        self.nautical_directions = {}
        for i,n in enumerate(self.names):
            self.paths[n] = functions.path_era([self.data['ra'][i], self.data['dec'][i]], self.observatory, self.era)
            self.lunar_distances[n] = functions.angle_to_moon([self.data['ra'][i], self.data['dec'][i]], self.lunar_coords)
            self.nautical_directions[n] = functions.calculate_pointing([self.data['ra'][i], self.data['dec'][i]], self.observatory, self.era)

    #Plotting function
    def plot_paths(self, targets = None, priorities = None, moon = True, lunar_distance = False, angle_limits = True, figsize = (15,10), title = None, angle_ax = False, dashapp = False, dark = False):
        fig = go.Figure()

        if dark:
            #fig.update_layout(plot_bgcolor="#444444", paper_bgcolor="#444444")
            fig.update_layout(plot_bgcolor="#1e1e1e", paper_bgcolor="#1e1e1e")
            colours = {'moon':'#ffffff', 'twilight':'rgba(255, 255, 255, 0.2)'}
        else:
            colours = {'moon':'#000000', 'twilight':'rgba(0, 0, 0, 0.2)'}
        if dashapp:
            fig.update_layout(height = 600)
    
        #Plotting specified target traces
        if targets != None:
            for j, name in enumerate(targets):
                fig.add_trace(go.Scatter(x=self.utc, y=self.paths[name], mode = 'lines', name = name))

                if lunar_distance:
                    for i in range(len(self.era)):
                        if self.twilights[2][0]-10 < self.era[i] < self.twilights[2][1]+10 and (i + 4*j) % 30 == 0 and 5 < self.paths[name][i] < 85:
                            fig.add_annotation(x=self.utc[i], y=self.paths[name][i], text=str(int(round(self.lunar_distances[name][i], 0))), showarrow=False, yshift=0)
        #Plotting traces of targets with specified priorities
        elif priorities != None:
            j = 0
            for name in self.data.index:
                if self.data.priority[name] in priorities:
                    fig.add_trace(go.Scatter(x=self.utc, y=self.paths[name], mode = 'lines', name = name))

                    if lunar_distance:
                        for i in range(len(self.era)):
                            if self.twilights[2][0]-10 < self.era[i] < self.twilights[2][1]+10 and (i + 4*j) % 30 == 0 and 5 < self.paths[name][i] < 85:
                                fig.add_annotation(x=self.utc[i], y=self.paths[name][i], text=str(int(round(self.lunar_distances[name][i], 0))), showarrow=False, yshift=0)
                    j += 1
        #Plotting all target traces
        else:
            for j, name in enumerate(self.data.index):
                fig.add_trace(go.Scatter(x=self.utc, y=self.paths[name], mode = 'lines', name = name))

                if lunar_distance:
                    for i in range(len(self.era)):
                        if self.twilights[2][0]-10 < self.era[i] < self.twilights[2][1]+10 and (i + 4*j) % 30 == 0 and 5 < self.paths[name][i] < 85:
                            fig.add_annotation(x=self.utc[i], y=self.paths[name][i], text=str(int(round(self.lunar_distances[name][i], 0))), showarrow=False, yshift=0)
            
        #Plotting moon trace
        if moon:
            fig.add_trace(go.Scatter(x=self.utc, y=self.lunar_path, mode = 'lines', name = 'Moon', line = dict(color = colours['moon'], dash = 'dot')))
        
        #Setting axis ranges
        fig.update_layout(yaxis_range = [-10, 90], yaxis_title = 'Elevation angle')
        fig.update_layout(xaxis_range = [functions.era_to_utc(self.twilights[2][0]-10), functions.era_to_utc(self.twilights[2][1]+10)], xaxis_title = 'UT')
    
        #Setting plot font
        if dashapp:
            if dark:
                fig.update_layout(font = dict(family="Trebuchet MS", size=16, color = '#dddddd'))
                fig.update_xaxes(gridcolor = '#666666')
                fig.update_yaxes(gridcolor = '#666666')
            else:
                fig.update_layout(font = dict(family="Trebuchet MS", size=16))
        else:   
            fig.update_layout(font = dict(family="Times",size=16))
    
        #Shading limiting angles
        if angle_limits:
            fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[self.min_angle, self.min_angle], line = dict(color = 'rgba(214, 27, 27, 0)'), fillcolor = 'rgba(214, 27, 27, 0.3)', showlegend=False, fill = 'tozeroy'))
            fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[100, 100], line = dict(color = 'rgba(214, 27, 27, 0)'), fillcolor = 'rgba(214, 27, 27, 0.3)', showlegend=False))
            fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[self.max_angle, self.max_angle], line = dict(color = 'rgba(214, 27, 27, 0)'), fillcolor = 'rgba(214, 27, 27, 0.3)', showlegend=False, fill = 'tonexty'))
        
        #Shading ground
        fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[0, 0], line = dict(color = 'black'), showlegend=False))
        fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[-90, -90], line = dict(color = 'black'), fill='tonexty', fillcolor='rgba(0, 0, 0, 0.8)', showlegend=False))
        
        #Shading sunsets
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[0][0])], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[0][0])], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Astronomical sunset'))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[1][0])], y=[100, 100], line = dict(color = 'black'), showlegend=False, marker={'opacity': 0}))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[1][0])], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Nautical sunset'))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[2][0])], y=[100, 100], line = dict(color = 'black'), showlegend=False, marker={'opacity': 0}))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[2][0])], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Civil sunset'))
        
        #Shading sunrises
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[0][1]), self.utc[-1]], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[0][1]), self.utc[-1]], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Astronomical sunrise'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[1][1]), self.utc[-1]], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[1][1]), self.utc[-1]], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Nautical sunrise'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[2][1]), self.utc[-1]], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[2][1]), self.utc[-1]], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Civil sunrise'))
        
        if dashapp:
            return fig
        else:
            fig.show()

    #Plot the derived queue
    def plot_queue(self, moon = True, angle_limits = True, figsize = (15,10), title = None, angle_ax = False, priority_colours = {1:'#1ABC9C', 2:'#FFB900', 3:'#FF9600', 4:'#D62727'}, dashapp = False, dark = False):
        fig = go.Figure()
        
        if dark:
            #fig.update_layout(plot_bgcolor="#444444", paper_bgcolor="#444444")
            fig.update_layout(plot_bgcolor="#1e1e1e", paper_bgcolor="#1e1e1e")
            colours = {'moon':'#ffffff', 'twilight':'rgba(255, 255, 255, 0.2)', 'path':'rgba(255, 255, 255, 0.1)'}
        else:
            colours = {'moon':'#000000', 'twilight':'rgba(0, 0, 0, 0.2)', 'path':'rgba(0, 0, 0, 0.1)'}
        if dashapp:
            fig.update_layout(height = 600)

        #Cutting path slices for the observations
        cut_paths = {}
        for name in self.queue.index:
            if '__gap' not in name:
                cut_paths[name] = functions.cut_array(self.era, self.paths[name], [self.queue['start_angle'][name], self.queue['start_angle'][name]+self.queue['obs_time_angles'][name]], z = self.utc)
                fig.add_trace(go.Scatter(x=self.utc, y=self.paths[name], mode = 'lines', name = name, line = dict(color = colours['path']), showlegend = False))
    
        #Plotting target traces
        for name in self.queue.index:
            if '__gap' not in name:
                fig.add_trace(go.Scatter(x=cut_paths[name][2], y=cut_paths[name][1], mode = 'lines', name = name, fill='tozeroy', line = dict(color = priority_colours[self.queue.priority[name]])))
    
        #Plotting moon trace
        if moon:
            fig.add_trace(go.Scatter(x=self.utc, y=self.lunar_path, mode = 'lines', name = 'Moon', line = dict(color = colours['moon'], dash = 'dot')))
        
        #Setting axis ranges
        fig.update_layout(yaxis_range = [-10, 90], yaxis_title = 'Elevation angle')
        fig.update_layout(xaxis_range = [functions.era_to_utc(self.twilights[2][0]-10), functions.era_to_utc(self.twilights[2][1]+10)], xaxis_title = 'UT')
    
        #Setting plot font
        if dashapp:
            if dark:
                fig.update_layout(font = dict(family="Trebuchet MS", size=16, color = '#dddddd'))
                fig.update_xaxes(gridcolor = '#666666')
                fig.update_yaxes(gridcolor = '#666666')
            else:
                fig.update_layout(font = dict(family="Trebuchet MS", size=16))
        else:   
            fig.update_layout(font = dict(family="Times",size=16))
        
        #Shading limiting angles
        if angle_limits:
            fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[self.min_angle, self.min_angle], line = dict(color = 'rgba(214, 27, 27, 0)'), fillcolor = 'rgba(214, 27, 27, 0.3)', showlegend=False, fill = 'tozeroy'))
            fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[100, 100], line = dict(color = 'rgba(214, 27, 27, 0)'), fillcolor = 'rgba(214, 27, 27, 0.3)', showlegend=False))
            fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[self.max_angle, self.max_angle], line = dict(color = 'rgba(214, 27, 27, 0)'), fillcolor = 'rgba(214, 27, 27, 0.3)', showlegend=False, fill = 'tonexty'))
    
        #Shading ground
        fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[0, 0], line = dict(color = 'black'), showlegend=False))
        fig.add_trace(go.Scatter(x=[self.utc[0], self.utc[-1]], y=[-90, -90], line = dict(color = 'black'), fill='tonexty', fillcolor='rgba(0, 0, 0, 0.8)', showlegend=False))
        
        #Shading sunsets
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[0][0])], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[0][0])], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Astronomical sunset'))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[1][0])], y=[100, 100], line = dict(color = 'black'), showlegend=False, marker={'opacity': 0}))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[1][0])], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Nautical sunset'))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[2][0])], y=[100, 100], line = dict(color = 'black'), showlegend=False, marker={'opacity': 0}))
        fig.add_trace(go.Scatter(x=[self.utc[0], functions.era_to_utc(self.twilights[2][0])], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Civil sunset'))
        
        #Shading sunrises
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[0][1]), self.utc[-1]], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[0][1]), self.utc[-1]], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Astronomical sunrise'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[1][1]), self.utc[-1]], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[1][1]), self.utc[-1]], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Nautical sunrise'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[2][1]), self.utc[-1]], y=[100, 100], line = dict(color = 'black'), showlegend=False, mode = 'lines'))
        fig.add_trace(go.Scatter(x=[functions.era_to_utc(self.twilights[2][1]), self.utc[-1]], y=[0, 0], line = dict(color = 'black'), fill='tonexty', fillcolor=colours['twilight'], showlegend=False, mode = 'lines', name = 'Civil sunrise'))
        
        if dashapp:
            return fig
        else:
            fig.show()

    def generate_queue(self, limiting_twilight = 'astronomical', restricted_directions = []):
        try:
            self.restricted_directions = []
            for d in restricted_directions:
                self.restricted_directions.append(d.upper())
    
            #Deciding which sunset to start the queue from
            if limiting_twilight == 'astronomical':
                self.limiting_twilight = self.twilights[0]
            elif limiting_twilight == 'nautical':
                self.limiting_twilight = self.twilights[1]
            elif limiting_twilight == 'civil':
                self.limiting_twilight = self.twilights[2]
            else:
                print(f"'{limiting_twilight}' not and acceptable twilight setting. The options are 'astronomical', 'nautical', and 'civil'. Defaulting to 'astronomical'.")
                self.limiting_twilight = self.twilights[0]
    
            #Creating empty dataframe to hold to queue
            self.queue = pd.DataFrame(columns = ['name', 'ra', 'dec', 'priority', 'obs_times', 'obs_time_angles', 'night_peak', 'start_angle', 'start_ut']).set_index('name')
            self.reserve = []
    
            #Setting data dataframe as variable within function accessible by duckdb SQL queries
            data = self.data_duckdb
            
            #Dividing up the target table into priorities ordered by peak ra throughout the night
            self.priorities = set(self.data['priority'])
            self.priority_tables = {}
            for p in self.priorities:
                self.priority_tables[p] = duckdb.query(f"SELECT * FROM data WHERE priority = {p} ORDER BY night_peak ASC;").df()
    
            for priority in self.priorities:
                for name in self.priority_tables[priority].name:
                    self.add_to_queue(name)
    
            #If we are yet to cross the morning twilight however have now observable targets remaining, introduce a gap to the queue
            try:
                while self.queue.start_angle[-1] + self.queue.obs_time_angles[-1] < self.limiting_twilight[1]:
                    self.add_gap_to_queue()
                    self.check_reserve_list()
            except:
                self.add_gap_to_queue()
                self.check_reserve_list()
                while self.queue.start_angle[-1] + self.queue.obs_time_angles[-1] < self.limiting_twilight[1]:
                    self.add_gap_to_queue()
                    self.check_reserve_list()

        except:
            pass 

    def add_to_queue(self, name):
        #Check if object meets the angular requirements somewhere between astronomical sunset and sunrise
        if self.check_valid_object(name):

            #If this is the first object add it at astronomical sunset
            if len(self.queue) == 0:
                start1 = int(np.ceil(self.limiting_twilight[0]))

                #Compiling row dictionary adding target after
                new_row1 = {'ra':self.data.ra[name], 'dec':self.data.dec[name], 'priority':self.data.priority[name], 'obs_times':self.data.obs_times[name], 'obs_time_angles':self.data.obs_time_angles[name], 'night_peak':self.data.night_peak[name], 'start_angle':start1, 'start_ut':functions.era_to_utc(start1)} 

            #Else calculate two start options
            #1 appending the target to the end
            #2 appending the target before the previously final target
            else:
                start1 = self.queue.start_angle[-1] + self.queue.obs_time_angles[-1]
                start2 = self.queue.start_angle[-1]
                replace_start2 = self.queue.start_angle[-1] + self.data['obs_time_angles'][name]
            
                #Compiling row dictionary adding target after
                new_row1 = {'ra':self.data.ra[name], 'dec':self.data.dec[name], 'priority':self.data.priority[name], 'obs_times':self.data.obs_times[name], 'obs_time_angles':self.data.obs_time_angles[name], 'night_peak':self.data.night_peak[name], 'start_angle':start1, 'start_ut':functions.era_to_utc(start1)} 
                new_row2 = {'ra':self.data.ra[name], 'dec':self.data.dec[name], 'priority':self.data.priority[name], 'obs_times':self.data.obs_times[name], 'obs_time_angles':self.data.obs_time_angles[name], 'night_peak':self.data.night_peak[name], 'start_angle':start2, 'start_ut':functions.era_to_utc(start2)} 
    
            #Creating new queue instance with new object appended at the end
            new_queue1 = self.queue.copy()
            new_queue1.loc[name] = new_row1
    
            if len(self.queue) == 0:
                #Checking valid queue
                valid1 = self.check_valid_queue(new_queue1)
                if valid1:
                    self.queue = new_queue1
                    if name in self.reserve:
                        self.reserve.remove(name)
                    self.check_reserve_list()
                    return
                else:
                    if name not in self.reserve:
                        self.reserve.append(name)
                    return
            
            else:
                #Creating new queue instance with new target appended before the previously final target
                new_queue2 = self.queue.copy()
                new_queue2['start_angle'][-1] = replace_start2
                new_queue2.loc[name] = new_row2
                
                #Checking the validity of the two queues
                valid1 = self.check_valid_queue(new_queue1)
                valid2 = self.check_valid_queue(new_queue2)
    
                if valid1 == True and valid2 == True:
                    #Both cases are valid, chose the queue that maximises the area underneath
                    new_queue_max = self.calculate_maximum_areas(new_queue1, new_queue2, [self.queue.index[-1], name])

                    self.queue = new_queue_max
                    if name in self.reserve:
                        self.reserve.remove(name)
                    self.check_reserve_list()
                    return
                elif valid1 == True and valid2 == False:
                    #Only new_queue1 is valid
                    self.queue = new_queue1
                    if name in self.reserve:
                        self.reserve.remove(name)
                    self.check_reserve_list()
                    return
                elif valid1 == False and valid2 == True:
                    #Only new_queue2 is valid
                    self.queue = new_queue2.sort_values('start_angle')
                    if name in self.reserve:
                        self.reserve.remove(name)
                    self.check_reserve_list()
                    return
                else:
                    if name not in self.reserve:
                        self.reserve.append(name)
                    return

    #Passing through the reserve list trying to add objects
    def check_reserve_list(self):
        #Reserve is empty, skip
        if len(self.reserve) == 0:
            return
        #Pass through reserves and attempt to add the to the queue
        else:
            for name in self.reserve:
                self.add_to_queue(name)

    #Running through the queue to test if the objects are observable
    def check_valid_queue(self, queue):
        for name in queue.index:
            if '__gap' not in name:
                cut_path = functions.cut_array(self.era, self.paths[name], [queue['start_angle'][name], queue['start_angle'][name]+queue['obs_time_angles'][name]], self.lunar_distances[name])
                cut_nautical_direction = functions.cut_array(self.era, self.nautical_directions[name], [queue['start_angle'][name], queue['start_angle'][name]+queue['obs_time_angles'][name]])[1]

                for i in range(len(cut_path[0])):
                    if self.max_angle < cut_path[1][i] or self.min_angle > cut_path[1][i] or cut_path[2][i] < self.minimum_lunar_distance or cut_nautical_direction[i] in self.restricted_directions:
                        #Invalid angle instance, lunar distance, or restricted pointing located, return False
                        return False
                
                    #Overlapping with final twilight
                    if cut_path[0][i] > self.limiting_twilight[1]:
                        return False

        #If no objects have thrown up a false flag, return True
        return True

    def check_valid_object(self, name):
        #Cut object path to between astronomical sunset and sunrise
        cut_path = functions.cut_array(self.era, self.paths[name], [self.limiting_twilight[0], self.limiting_twilight[1]], self.lunar_distances[name])
        cut_nautical_direction = functions.cut_array(self.era, self.nautical_directions[name], [self.limiting_twilight[0], self.limiting_twilight[1]])[1]

        for i in range(len(cut_path[0])):
            if self.min_angle < cut_path[1][i] < self.max_angle and cut_path[2][i] > self.minimum_lunar_distance and cut_nautical_direction[i] not in self.restricted_directions:
                #Valid angle instance, lunar distance, and poiting located, return True
                return True
        return False

    #Comparing the total areas from two queues and returning the maximum
    def calculate_maximum_areas(self, queue1, queue2, objs):
        #If one of the two objects is a gap, skip the check
        for name in objs:
            if '__gap' in name:
                return queue1

        cut_paths1 = []
        cut_paths2 = []

        #Cutting the paths according to the two queues
        for name in objs:
            cut_paths1.append(functions.cut_array(self.era, self.paths[name], [queue1['start_angle'][name], queue1['start_angle'][name]+queue1['obs_time_angles'][name]]))
            cut_paths2.append(functions.cut_array(self.era, self.paths[name], [queue2['start_angle'][name], queue2['start_angle'][name]+queue2['obs_time_angles'][name]]))

        #Initialising the queue areas
        area1 = 0
        area2 = 0

        #Calculating the area under the observed curves
        for i in range(len(cut_paths1)):
            area1 += np.trapz(cut_paths1[i][1], cut_paths1[i][0])
            area2 += np.trapz(cut_paths2[i][1], cut_paths2[i][0])

        #If plaing the new target before the previous target, return the corresponding queue
        if area2 > area1:
            return queue2.sort_values('start_angle')
        #Otherwise return the original queue
        else:
            return queue1

    #Function to introduce a gap (or extend an existing gap) in a queue due to a lack of observable objects
    def add_gap_to_queue(self):
        try:
            #If the most recent item in the queue is a gap, extend the gap
            if '__gap' in self.queue.index[-1]:
                self.queue['obs_time_angles'][-1] += 1
            #Otherwise introduce a new gap
            else:
                index = 1
                for name in self.queue.index:
                    if '__gap' in name:
                        index += 1
                new_row = {'ra':None, 'dec':None, 'priority':None, 'obs_times':None, 'obs_time_angles':1, 'night_peak':None, 'start_angle':self.queue['start_angle'][-1] + self.queue['obs_time_angles'][-1], 'start_ut':functions.era_to_utc(self.queue['start_angle'][-1] + self.queue['obs_time_angles'][-1])} 
                self.queue.loc[f'__gap{index}__'] = new_row
        except:
            index = 1
            new_row = {'ra':None, 'dec':None, 'priority':None, 'obs_times':None, 'obs_time_angles':1, 'night_peak':None, 'start_angle':self.limiting_twilight[0], 'start_ut':functions.era_to_utc(self.limiting_twilight[0])} 
            self.queue.loc[f'__gap{index}__'] = new_row
