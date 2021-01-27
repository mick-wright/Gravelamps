TARGET = gwlensing
OBJECTS = plCalc
LENSING_FOLDER = gwlensing/lensing
INFERENCE_FOLDER = gwlensing/inference

install: $(OBJECTS)
	@echo "Installing $(TARGET)"
	@pip install .
	@echo "$(TARGET) installed!" 

uninstall: clean_plCalc
	@echo "Uninstalling $(TARGET)"
	@pip uninstall $(TARGET)
	@echo "Uninstallation complete" 

plCalc:
	@echo "Building plCalc"
	$(MAKE) -C $(LENSING_FOLDER) copy
	@echo "plCalc built!" 

clean_plCalc:
	@echo "Removing plCalc"
	$(MAKE) -C $(LENSING_FOLDER) clean
	@echo "plCalc removed" 
