import plotly.express as px
import plotly.graph_objects as go

import numpy as np

def get_slider_no_steps(transition_dur=0, prefix=None):
    slider_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "visible": False,
            "xanchor": "right"
        },
        "transition": {"duration": transition_dur, "easing":"quad-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }
    if prefix:
        slider_dict["currentvalue"]["visible"] = True
        slider_dict["currentvalue"]["prefix"] = prefix
    return slider_dict.copy()

def get_slider_step(frame_nr, frame_duration=500):
    slider_step = {"args": [
        [str(frame_nr)],
        {"frame": {"duration": frame_duration, "redraw": True},
         "mode": "immediate",
         "transition": {"duration": 0}}
         ],
        "label": str(frame_nr),
        "method": "animate"}
    return slider_step

def add_animation_frame(data, frame_nr):
    fr = go.Frame(data=data, name=str(frame_nr))
    return fr

def get_pause_play(frame_duration=500, transition_dur = 0):
    update_menues = {
        "buttons": [
            {
                "args": [None, {"frame": {"duration": frame_duration, "redraw": True},
                                "fromcurrent": True, "transition": {"duration": transition_dur, "easing":"quad-in-out"}}],
                "label": "Play",
                "method": "animate"
            },
            {
                "args": [[None], {"frame": {"duration": 10, "redraw": True},
                                  "mode": "immediate",
                                  "transition": {"duration": transition_dur}}],
                "label": "Pause",
                "method": "animate"
            }
        ],
        "direction": "left",
        "pad": {"r": 10, "t": 87},
        "showactive": False,
        "type": "buttons",
        "x": 0.1,
        "xanchor": "right",
        "y": 0,
        "yanchor": "top"
    }
    return update_menues

