def levels(SiteID):
    """ Dictionary of all of the sites with different levels other than 0 at the surface.
	Arg: 
	SiteID (str) : The Site ID you input to retrieve the levels.
    
    A.F. now adapted for WACCM 70 L
    
    """
    level = {
    'ZEP': 2,
    'AND': 1,
    'MWA': 3,
    'MLO': 9,
    'LLN':8,
        }
    return level.get(SiteID.upper(), 0)
