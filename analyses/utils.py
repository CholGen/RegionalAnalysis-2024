import matplotlib.dates as mdates
from matplotlib.collections import LineCollection
import matplotlib as mpl
from datetime import datetime as dt
from datetime import timedelta
import time
from subprocess import run

import pandas as pd
import numpy as np


def get_okabe_ito_palette():
    return ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


def setup_plotting_standards():
    prop = mpl.font_manager.FontProperties( 'Roboto' )
    mpl.rcParams['font.sans-serif'] = prop.get_name()
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.weight'] = 300
    mpl.rcParams['axes.labelweight'] = 300
    mpl.rcParams['font.size'] = 16

    mpl.rcParams['figure.dpi'] = 200

    COLOR = '#343434'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR


def timeseries_formatting( ax ):
    """
    Specifies the proper date formatting of the xaxis.
    """
    # Properly label timeseries
    ax.xaxis.set_minor_locator( mdates.MonthLocator() )
    ax.xaxis.set_minor_formatter( mdates.DateFormatter( '%b' ) )
    ax.xaxis.set_major_locator( mdates.YearLocator() )
    ax.xaxis.set_major_formatter( mdates.DateFormatter( '%b\n%Y' ) )


def skipped_timeseries_formatting( ax ):
    ax.xaxis.set_major_locator( mdates.MonthLocator() )
    formatter = mdates.DateFormatter( '%b' )
    long_formatter = mdates.DateFormatter( '%b\n%Y' )
    accepted_months = ["Jan", "Mar", "May", "Jul", "Sep", "Nov"]
    labels = []
    for date in ax.get_xticks():
        formatted_date = formatter.format_data( date )
        if formatted_date in accepted_months:
            if formatted_date == "Jan":
                labels.append( long_formatter.format_data( date ) )
            else:
                labels.append( formatted_date )
        else:
            labels.append( "" )
    ax.set_xticklabels( labels )


def basic_formatting( ax, spines=["left", "bottom"], which="y", title=None, ylabel=None, xlabel=None, xsize=12, ysize=12,
                      xlims=None, ylims=None ):
    """
    Specifies the basic formatting of all plots.
    """
    # Remove spines
    [ax.spines[j].set_visible( False ) for j in ax.spines if j not in spines]

    # Format axes ticks
    ax.tick_params( axis="x", which="both", direction="inout", labelbottom=True, labelsize=xsize )
    ax.tick_params( axis="y", which="both", direction="inout", labelleft=True, labelsize=ysize )

    # Label axes
    ax.set_xlabel( xlabel, fontsize=xsize )
    ax.set_ylabel( ylabel, fontsize=ysize )
    ax.set_title( title )

    ax.set_facecolor( "w" )

    # Add a simple grid
    ax.grid( which="both", axis=which, linewidth=1, color="#F1F1F1", zorder=1 )

    # Add the xlims
    if xlims:
        ax.set_xlim( xlims )
    if ylims:
        ax.set_ylim( ylims )


def _toYearFraction( date, date_format="%Y-%m-%d" ):
    """ Converts datetime object to a decimal year
    Parameters
    ----------
    date: datetime.datetime
        date to be converted.

    Returns
    -------
    float
        date in decimal year format.
    """

    def sinceEpoch( d ):  # returns seconds since epoch
        return time.mktime( d.timetuple() )

    date = dt.strptime( date, date_format )

    year = date.year
    start_of_this_year = dt( year=year, month=1, day=1 )
    start_of_next_year = dt( year=year + 1, month=1, day=1 )

    year_elapsed = sinceEpoch( date ) - sinceEpoch( start_of_this_year )
    year_duration = sinceEpoch( start_of_next_year ) - sinceEpoch( start_of_this_year )
    fraction = year_elapsed / year_duration

    return date.year + fraction


def dec_to_date( date_str ):
    """
    Converts a decimal date into a datetime object.
    Useful for interacting with output from BEAST.
    """
    if date_str is None:
        return None
    year = int( date_str )
    rem = date_str - year
    base = dt( year, 1, 1 )
    result = base + timedelta( seconds=(base.replace( year=base.year + 1 ) - base).total_seconds() * rem )
    return result

def hpd( data, level ):
    """
    Return highest posterior density interval from a list,
    given the percent posterior density interval required.
    """
    d = list( data )
    d.sort()

    nData = len( data )
    nIn = int( round( level * nData ) )
    if nIn < 2:
        return None
    # raise RuntimeError("Not enough data. N data: %s"%(len(data)))

    i = 0
    try:
        r = d[i + nIn - 1] - d[i]
    except IndexError:
        print( i )
        print( nIn )
        print( d )
        raise

    for k in range( len( d ) - (nIn - 1) ):
        rk = d[k + nIn - 1] - d[k]
        if rk < r:
            r = rk
            i = k

    assert 0 <= i <= i + nIn - 1 < len( d )

    return (d[i], d[i + nIn - 1])

def extract_metadata(
        alignment = "../data/combined_alignment.bcf.gz",
        cholgen_md = "../data/combined_cholgen_metadata.csv.gz",
        background_md = '../data/vc_metadata_workshop.csv',
        columns=None,
        rename=None ) -> pd.DataFrame:


    if columns is None:
        columns = ["terraID", "collection_date", "geo_loc_country"]
    if rename is None:
        rename = ["taxa", "collection_year", "country"]

    query = run( f"bcftools query -l {alignment}", shell=True, capture_output=True, text=True )
    assert query.returncode == 0, "BCFTools query was not successful, check output."
    taxa = query.stdout.strip().split( "\n" )

    cg = pd.read_csv( cholgen_md, compression="gzip" )
    cg = cg[columns]
    cg.columns = rename
    cg = cg.loc[~cg["taxa"].isna()]
    cg["origin"] = "CholGen"

    md = pd.read_csv( background_md )
    md["origin"] = "background"
    md = pd.concat( [md, cg] )
    md["continent"] = md["continent"].fillna( "Africa" )
    md["te"] = md["te"].fillna( "UND" )
    md["ref"] = md["ref"].fillna( "?" )
    md["workshop"] = md["workshop"].fillna( True )
    md = md.loc[md["taxa"].isin( taxa )].copy()
    md = md.drop_duplicates( subset=["taxa"], keep="first" )
    md["country"] = md["country"].replace( {"DRC": "Democratic Republic of the Congo"} )
    md["collection_year"] = pd.to_datetime( md["collection_year"], format="ISO8601" )
    md["year"] = md["collection_year"].dt.year

    return md

def load_africa_shapefile():
    afr = gpd.read_file( "/Users/natem/Documents/Data/shapefiles/Africa_admin0/afr_g2014_2013_0.shp" )
    afr = afr[["ADM0_NAME", "ISO3", "geometry"]]
    afr = afr.loc[
        ~afr["ADM0_NAME"].isin( ["Seychelles", "Sao Tome and Principe", "Cape Verde", "Comoros", "Mauritius"] )]
    afr = afr.explode()
    afr["area"] = afr["geometry"].area
    afr = afr.sort_values( "area", ascending=False ).groupby( "ADM0_NAME" ).first()
    afr = afr.reset_index()
    afr["geometry"] = afr["geometry"].simplify( 0.1 )
    return afr

# This is a rip of @evogytis baltic plotting function, but with the ability to plot horizontally.
def plotTree( self, ax, connection_type=None, target=None,
              x_attr=None, y_attr=None, width=None,
              colour=None, horizontal=False, **kwargs ):
    if target == None: target = lambda k: True
    if x_attr == None: x_attr = lambda k: k.x
    if y_attr == None: y_attr = lambda k: k.y
    if width == None: width = 2
    if colour == None: colour = 'k'
    if connection_type == None: connection_type = 'baltic'
    assert connection_type in ['baltic', 'direct', 'elbow'], 'Unrecognised drawing type "%s"' % (tree_type)

    branches = []
    colours = []
    linewidths = []
    for k in filter( target, self.Objects ):  ## iterate over branches
        x = x_attr( k )  ## get branch x position
        xp = x_attr( k.parent ) if k.parent else x  ## get parent x position
        y = y_attr( k )  ## get y position

        try:
            colours.append( colour( k ) ) if callable( colour ) else colours.append( colour )
        except KeyError:
            colours.append( (0.7, 0.7, 0.7) )
        linewidths.append( width( k ) ) if callable( width ) else linewidths.append( width )

        if connection_type == 'baltic':
            if horizontal:
                branches.append( ((y, xp), (y, x)) )
            else:
                branches.append( ((xp, y), (x, y)) )

            if k.is_node():
                yl, yr = y_attr( k.children[0] ), y_attr( k.children[-1] )
                if horizontal:
                    branches.append( ((yl, x), (yr, x)) )
                else:
                    branches.append( ((x, yl), (x, yr)) )
                linewidths.append( linewidths[-1] )
                colours.append( colours[-1] )
        elif connection_type == 'elbow':
            yp = y_attr( k.parent ) if k.parent else y  ## get parent x position
            branches.append( ((xp, yp), (xp, y), (x, y)) )
        elif connection_type == 'direct':
            yp = y_attr( k.parent )  ## get y position
            branches.append( ((xp, yp), (x, y)) )
        else:
            pass  ## for now

    if 'capstyle' not in kwargs: kwargs['capstyle'] = 'projecting'
    line_segments = LineCollection( branches, lw=linewidths, color=colours, **kwargs )
    ax.add_collection( line_segments )
    return ax

# This is a rip of @evogytis baltic plotting function, but with the ability to plot horizontally.
def plotPoints( tree, ax, x_attr=None, y_attr=None, target=None, size=None, colour=None,
                zorder=None, outline=None, outline_size=None, outline_colour=None, horizontal=False, style="baltic", linewidth=1, **kwargs ):
    if target == None: target = lambda k: k.is_leaf()
    if x_attr == None: x_attr = lambda k: k.x
    if y_attr == None: y_attr = lambda k: k.y
    if size == None: size = 40
    if colour == None: colour = lambda f: 'k'
    if zorder == None: zorder = 3

    if outline == None: outline = True
    if outline_size == None: outline_size = lambda k: size( k ) * 2 if callable( size ) else size * 2
    if outline_colour == None: outline_colour = 'k'

    xs = []
    ys = []
    colours = []
    sizes = []

    outline_xs = []
    outline_ys = []
    outline_colours = []
    outline_sizes = []
    for k in filter( target, tree.Objects ):

        xs.append( x_attr( k ) )
        ys.append( y_attr( k ) )

        colours.append( colour( k ) ) if callable( colour ) else colours.append( colour )
        sizes.append( size( k ) ) if callable( size ) else sizes.append( size )

        if outline:
            outline_xs.append( xs[-1] )
            outline_ys.append( ys[-1] )
            outline_colours.append( outline_colour( k ) ) if callable( outline_colour ) else outline_colours.append(
                outline_colour )
            outline_sizes.append( outline_size( k ) ) if callable( outline_size ) else outline_sizes.append(
                outline_size )

    if horizontal:
        xs, ys = ys, xs
        outline_xs, outline_ys = outline_ys, outline_xs

    if style == "baltic":
        ax.scatter( xs, ys, s=sizes, facecolor=colours, edgecolor='none', zorder=zorder,
                    **kwargs )  ## put a circle at each tip
        if outline:
            ax.scatter( outline_xs, outline_ys, s=outline_sizes, facecolor=outline_colours, edgecolor='none',
                        zorder=zorder - 1, **kwargs )  ## put a circle at each tip
    elif style == "nate":
        ax.scatter( xs, ys, s=sizes, facecolor=colours, edgecolor=outline_colours, linewidth=linewidth, zorder=zorder, **kwargs )

    return ax

def mm( x, V: float, M: float ):
	return ( V * x ) / ( M + x )

def weibull( x, V: float, b: float, c: float ) -> float:
    return V / (1 + np.exp( (b - np.log( x )) / c) )