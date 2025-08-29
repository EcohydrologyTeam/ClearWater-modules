from config import init_from_file

config = "./clearwater_modules_v2/dev.yml"

### Define the model
model = init_from_file(config)

model.run()

prt = 1
