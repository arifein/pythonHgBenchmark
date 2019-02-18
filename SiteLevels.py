def levels(SiteID):
    """ Dictionary of all of the sites with different levels other than 0 at the surface.
	Arg: 
	SiteID (str) : The Site ID you input to retrieve the levels.
    
    
    """
    level = {
    'ZEP': 3,
    'AND': 2,
    'MWA': 1,
    'MLO': 18,
    'MBO': 16,
    'NamCo': 18,
    'LLN':16,
        }
    return level.get(SiteID.upper(), 0)
