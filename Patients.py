import re


class Patient:
    """
    Used to represent information about each patient
    Gender, ancestry....
    """
    def __init__(self, line):
        """
        Patients should be read in line by line from a csv file
        formatted like this:
        PatientName,    Gender, Population, Superpopulation
        HG00096,        1,      GBR,        EUR
        Read each line into Patient()
        """
        self.patient_name, gender, self.population, self.super_population = line.split(",")
        if gender == "1":
            self.gender = "M"
        elif gender == "2":
            self.gender = "F"

    def __hash__(self):
        return hash(self.patient_name)

    def __eq__(self, other):
        return self.patient_name == other.patient_name


class Patients:
    def __init__(self, subject_info_file_path):
        """
        Creates a list of Patients from an input file
        """
        patient_match = re.compile("[A-Z]{2}\d{5}")
        self.patients = {}
        with open(subject_info_file_path) as subject_info:
            for line in subject_info:
                if patient_match.match(line) is not None:
                    p = Patient(line.strip())
                    self.patients[p.patient_name] = p

    def __getitem__(self, item):
        return self.patients[item]

    def __contains__(self, item):
        return item in self.patients
