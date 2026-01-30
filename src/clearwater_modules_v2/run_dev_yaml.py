from config.init import init_from_file
import datetime

# config = r"C:\Users\ptomasula\Repositories\ClearWater-modules\src\clearwater_modules_v2\dev_nsm.yml"
config = r"C:\Users\ptomasula\Repositories\ClearWater-modules\src\clearwater_modules_v2\dev.yml"

### Define the model
start = datetime.datetime.now()
print("Building Model")
model = init_from_file(config)
print(f"Time to build model: {datetime.datetime.now() - start}")

print("Running Model")
start = datetime.datetime.now()
model.run()
print(f"Time to run: {datetime.datetime.now() - start}")

prt = 1
