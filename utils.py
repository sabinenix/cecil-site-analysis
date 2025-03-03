import cecil
import json
import geopandas as gpd
import pandas as pd
from shapely.geometry import shape



def show_linked_requests(dataset_mapping):
    """For the current organisation, list all the data requests with their 
    relevant reprojections, datasets, resolutions, and the AOI names and IDs."""
    client = cecil.Client()

    # List all reprojections, data requests and AOIs in an organisation.
    reprojections = client.list_reprojections()
    data_requests = client.list_data_requests()
    aois = client.list_aois()

    # Convert objects to dictionaries.
    reprojection_dicts = [vars(r) for r in reprojections]
    data_request_dicts = [vars(dr) for dr in data_requests]
    aoi_dicts = [vars(aoi) for aoi in aois]

    # Convert lists of dictionaries to DataFrames.
    reprojection_df = pd.DataFrame(reprojection_dicts)
    data_request_df = pd.DataFrame(data_request_dicts)
    aoi_df = pd.DataFrame(aoi_dicts)
    
    # Rename id columns before merging to clarify type of id.
    reprojection_df = reprojection_df.rename(columns={'id': 'reprojection_id'})
    data_request_df = data_request_df.rename(columns={'id': 'data_request_id'})
    aoi_df = aoi_df.rename(columns={'id': 'aoi_id', 'name': 'aoi_name'})

    # Add dataset names using the dataset_mapping dictionary.
    data_request_df['dataset_name'] = data_request_df['dataset_id'].map(dataset_mapping)
    
    # First merge data requests with AOIs.
    dr_aoi_df = data_request_df.merge(
        aoi_df,
        left_on='aoi_id',
        right_on='aoi_id',
        how='left'
    )
    
    # Next, merge the data requests with the reprojections.
    result_df = dr_aoi_df.merge(
        reprojection_df,
        left_on='data_request_id',
        right_on='data_request_id',
        how='left'  # Keep all data requests even if no reprojection exists
    )

    # Select and order relevant columns for resulting dataframe.
    result_df = result_df[['data_request_id', 'reprojection_id', 'dataset_name', 'resolution', 'aoi_name', 'aoi_id']]

    return result_df


def get_reprojection_details_by_id(reprojection_id, dataset_mapping):
    """For a given reprojection ID, return the dataset name, aoi name and ID."""
    
    # Get full dataframe linking reprojections to datasets and AOIs.
    result_df = show_linked_requests(dataset_mapping)

    # Filter to get just the reprojection ID.
    filtered = result_df[result_df['reprojection_id'] == reprojection_id]

    if filtered.empty:
        return {"error": "Reprojection ID not found"}

    # Extract details of interest for querying data.
    dataset_name = filtered['dataset_name'].iloc[0]
    aoi_name = filtered['aoi_name'].iloc[0]
    aoi_id = filtered['aoi_id'].iloc[0]

    return dataset_name, aoi_name, aoi_id

def get_data_request_details_by_id(data_request_id, dataset_mapping):
    """For a given data request ID, return the dataset name, aoi name and ID."""
    
    # Get full dataframe linking reprojections to datasets and AOIs.
    result_df = show_linked_requests(dataset_mapping)

    # Filter to get just the reprojection ID.
    filtered = result_df[result_df['data_request_id'] == data_request_id]

    if filtered.empty:
        return {"error": "Data Request ID not found"}

    # Extract details of interest for querying data.
    dataset_name = filtered['dataset_name'].iloc[0]
    aoi_name = filtered['aoi_name'].iloc[0]
    aoi_id = filtered['aoi_id'].iloc[0]

    return dataset_name, aoi_name, aoi_id


def query_all_data_analytics(reprojection_id, dataset_mapping):
    """Pull all data for a given reprojection_id."""

    client = cecil.Client()

    print(f'Attempting to query data for reprojection: {reprojection_id}')

    # Get details for querying the data from the reprojection_id.
    dataset_name, aoi_name, aoi_id = get_reprojection_details_by_id(reprojection_id, 
                                                          dataset_mapping)
    print(f'Pulling {dataset_name} for AOI {aoi_name} and AOI id {aoi_id}.')
    
    # Full table scan to pull all data for a given reprojection_id.
    df = client.query(f'''
        SELECT * 
        FROM 
            {dataset_name}
        WHERE
            aoi_id = '{aoi_id}' AND reprojection_id = '{reprojection_id}'
        ''')
    return df

def query_all_data_raw(data_request_id, dataset_mapping):
    """Pull all data from the raw_db for a given data_request_id."""

    client = cecil.Client()

    print(f'Attempting to query data for data request: {data_request_id}')

    # Get details for querying the data from the reprojection_id.
    dataset_name, aoi_name, aoi_id = get_data_request_details_by_id(data_request_id, 
                                                                    dataset_mapping)
    print(f'Pulling {dataset_name} for AOI {aoi_name} and AOI id {aoi_id}.')
    
    # Full table scan to pull all data for a given reprojection_id.
    df = client.query(f'''
        SELECT * 
        FROM 
            raw_db.{dataset_name}
        WHERE
            aoi_id = '{aoi_id}' 
            AND data_request_id = '{data_request_id}'
        ''')
    return df


def df_to_gdf(df, crs="EPSG:4326"):
    """
    Convert a DataFrame (from Snowflake) to GeoDataFrame with 'geometry' column.
    
    Args: 
    df: DataFrame pulled from Snowflake with 'pixel boundary' column.
    crs: (String) CRS of the data. Default is 'EPSG:4326'. 
    
    Returns:
    gdf: GeoDataFrame created using the 'geometry' column.
    """
    
    df["geometry"] = df["pixel_boundary"].apply(
    lambda x: shape(json.loads(x)) if isinstance(x, str) else None
    )
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=crs)
    return gdf
