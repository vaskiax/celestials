from pandas import DataFrame, Series
def OrbitalsQuery (astro:str, params:list[str]) -> DataFrame:
    """
    Función que obtiene los elementos orbitales de un cuerpo celeste.

    Parameters
    ----------
    astro : str
        Nombre del cuerpo celeste.
    params : list
        Lista con los parametros a obtener.

    Returns
    -------
    pd.DataFrame
        DataFrame con los elementos orbitales del cuerpo celeste.

    """
    from astroquery.jplsbdb import SBDB
    sbdb = DataFrame(SBDB.query(astro))

    return Series({f'{name}': sbdb.loc['elements', 'orbit'][name].value
                      if name != 'e' else sbdb.loc['elements', 'orbit'][name]
                       for name in params})


def MajorBodyQuery(id:int, params:list[str]) -> DataFrame:
    """
    Función que obtiene los elementos orbitales de un cuerpo celeste.

    Parameters
    ----------
    astro : str
        Nombre del cuerpo celeste.
    params : list
        Lista con los parametros a obtener.

    Returns
    -------
    pd.DataFrame
        DataFrame con los elementos orbitales del cuerpo celeste.

    """
    from astroquery.jplhorizons import Horizons
    obj = Horizons(id=id, location=None, epochs=None)
    obj_elements = obj.elements()

    return Series({f'{name}': obj_elements[name].value.data[0] for name in params})

