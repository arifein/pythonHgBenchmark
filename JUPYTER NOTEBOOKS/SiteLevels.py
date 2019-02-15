def levels(SiteID):
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
