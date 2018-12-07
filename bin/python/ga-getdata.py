from oauth2client.service_account import ServiceAccountCredentials
from apiclient.discovery import build
from datetime import datetime
from datetime import timedelta
import os
import sys

key_file = sys.argv[1]
if not os.path.isfile(key_file):
	print("Cannot find a key file at '" + key_file + "'")
	print("Please specify the path to key file as the first command line argument")
	sys.exit(1)

# The scope for the OAuth2 request.
SCOPE = 'https://www.googleapis.com/auth/analytics.readonly'

VIEW_ID = '51277087'

# Create service account credentials instance to authorize the GA API query
credentials = ServiceAccountCredentials.from_json_keyfile_name(key_file, SCOPE)

# Build the service object.
analytics = build('analyticsreporting', 'v4', credentials=credentials)

# Calculate the query date range
# isoWeen starts on Monday, so start date is on Monday
# and end date is on Sunday
endDate = datetime.now() - timedelta(days = datetime.now().weekday() + 1)
startDate = endDate - timedelta(weeks = 32) + timedelta(days=1)

# run the query
res = analytics.reports().batchGet(
      body={
        'reportRequests': [
        {
          'viewId': VIEW_ID,
          'dateRanges': [{'startDate': startDate.strftime('%Y-%m-%d'), 'endDate': endDate.strftime('%Y-%m-%d')}],
          'metrics': [{'expression': 'ga:totalEvents'}],
          'dimensions': [{'name': 'ga:isoYearIsoWeek'}],
          'dimensionFilterClauses': [{
	  	'filters' : [{
                      'dimensionName':'ga:eventAction',
                      'operator' : 'EXACT',
                      'expressions' : ['tide-search']}]
           }]
         },
         {
          'viewId': VIEW_ID,
          'dateRanges': [{'startDate': startDate.strftime('%Y-%m-%d'), 'endDate': endDate.strftime('%Y-%m-%d')}],
          'metrics': [{'expression': 'ga:totalEvents'}],
          'dimensions': [{'name': 'ga:isoYearIsoWeek'}]
         }
        ]
      }
  ).execute()

# parse the results

row_cnt = (len(res['reports'][0]['data']['rows']))
for i in range(0, row_cnt):
   week = res['reports'][0]['data']['rows'][i]['dimensions'][0]
   search_cnt = res['reports'][0]['data']['rows'][i]['metrics'][0]['values'][0]
   total_cnt = res['reports'][1]['data']['rows'][i]['metrics'][0]['values'][0]
   print(','.join([week, total_cnt, search_cnt]))

