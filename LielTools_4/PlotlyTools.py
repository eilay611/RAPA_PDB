
import plotly
import plotly.graph_objects as go
import plotly.express as px

def plot_treemap(df, path_vars, values_vars, color_vals, save_path):
    fig = px.treemap(df,
                     path=path_vars,
                     values=values_vars,
                     color=color_vals
                    )

    plotly.offline.plot(fig, filename=save_path)