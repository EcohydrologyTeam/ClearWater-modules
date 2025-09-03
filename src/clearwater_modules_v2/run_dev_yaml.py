from config import init_from_file

config = "./clearwater_modules_v2/dev.yml"

### Define the model
model = init_from_file(config)

import datetime

start = datetime.datetime.now()
model.run()

print(f"Time to run: {datetime.datetime.now() - start}")

prt = 1
