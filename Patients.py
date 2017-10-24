import re


class Patient:
    """
    Used to represent information about each patient
    Gender, ancestry....
    """
    def __init__(self, line=None, patient_name=None):
        """
        Two options;
        From csv file...
        Patients should be read in line by line from a csv file
        formatted like this:
        PatientName,    Gender, Population, Superpopulation
        HG00096,        1,      GBR,        EUR
        Read each line into Patient()

        From patient_name, leaving other fields blank
        """
        if line is not None:
            self.patient_name, gender, self.population, self.super_population = line.split(",")
            if gender == "1":
                self.gender = "M"
            elif gender == "2":
                self.gender = "F"
        if patient_name is not None:
            self.patient_name = patient_name
            self.gender = ""
            self.population = ""
            self.super_population = ""

    def __hash__(self):
        return hash(self.patient_name)

    def __eq__(self, other):
        return self.patient_name == other.patient_name


class Patients:
    def __init__(self, patients=None, subject_info_file_path=None):
        """
        Creates a list of Patients from an input file or from an existing list
        Call this only with patients OR subject_info_file_path
        """
        if patients is not None:
            self.patients = {}
            for patient_name in patients:
                self.patients[patient_name] = Patient(patient_name=patient_name)
        if subject_info_file_path is not None:
            self.patients = {}
            with open(subject_info_file_path) as subject_info:
                subject_info.readline()  # Get rid of header line
                for line in subject_info:
                    if bool(re.search(r'\d', line)):
                        p = Patient(line=line.strip())
                        self.patients[p.patient_name] = p

    def __getitem__(self, item):
        return self.patients[item]

    def __contains__(self, item):
        return item in self.patients
