"""WIP ROICaT processing task for MPCI data."""
import logging
from itertools import groupby
from collections import Counter

from scipy.cluster.hierarchy import linkage, fcluster
import numpy as np
from one.alf.path import ALFPath
import roicat
import roicat.data_importing

from ibllib.pipes.tasks import Task
from ibllib.oneibl.data_handlers import ExpectedDataset


logger = logging.getLogger('ibllib.ROICaT')


class SubjectAggregateTask(Task):
    """A base class for tasks that operate on a subject level, aggregating data across sessions."""

    def __init__(self, subject_path, **kwargs):
        """A task for running ROICaT accross a subject's imaging sessions.

        Parameters
        ----------
        subject_path : pathlib.Path
            The root subject path containing the sessions to process.
        """
        self.subject_path = subject_path
        super().__init__(self.subject_path, **kwargs)


class ROICaTTask(SubjectAggregateTask):
    cpu = 1   # CPU resource
    gpu = 1   # GPU resources: as of now, either 0 or 1
    io_charge = 30  # integer percentage
    priority = 70  # integer percentage, 100 means highest priority
    ram = 12  # RAM needed to run (GB)
    job_size = 'large'  # either 'small' or 'large', defines whether task should be run as part of the large or small job services
    env = 'roicat'  # the environment name within which to run the task (NB: the env is not activated automatically!)

    def __init__(self, subject_path, select_RFM_days_only=False, **kwargs):
        """A task for running ROICaT accross a subject's imaging sessions.

        Parameters
        ----------
        subject_path : pathlib.Path
            The root subject path containing the sessions to process.
        select_RFM_days_only : bool, optional
            Whether to only select days with Receptive Field Mapping data, by default False.
        """
        super().__init__(subject_path, **kwargs)
        self.select_RFM_days_only = select_RFM_days_only

    @property
    def signature(self):
        # The number of in and outputs will be dependent on the number of input raw imaging folders and output FOVs
        I = ExpectedDataset.input  # noqa
        # TODO Handle mlapdv provenance properly
        inputs = I('mpciMeanImage.mlapdv_estimate.npy', '????-??-??/???/alf/FOV_??', True, unique=False)
        inputs &= I('_suite2p_ROIData.raw.zip', '????-??-??/???/alf/FOV_??', True, unique=False)
        if self.select_RFM_days_only:
            raise NotImplementedError('AND with aggregates not implemented yet')
            # FIXME may be in 'alf/' if only one protocol run that session
            inputs &= I('_ibl_passiveRFM.times.npy', '????-??-??/???/alf/task_??', True, unique=False)
        # TODO define output files
        outputs = []
        return {'input_files': [inputs], 'output_files': outputs}

    def get_data_handler(self, **kwargs):
        """Data handlers not supported yet."""
        return None

    def _run(self, sessions_to_exclude=None):
        """Run ROICaT processing for the subject."""
        # Load the list of sessions to process
        paths = self.fetch_fov_list(sessions_to_exclude=sessions_to_exclude)
        clusters = self.group_fovs(paths)
        for cID, paths in clusters.items():
            # Load the data for ROICaT processing
            data = self.load_data(paths)
            
        
    def load_data(self, paths):
        """Load the data for ROICaT processing."""
        
        data = roicat.data_importing.Data_suite2p(
            paths_statFiles=...,
            paths_opsFiles=...,
            um_per_pixel=1.2,  ## IMPORTANT PARAMETER
            new_or_old_suite2p='new',
            type_meanImg='meanImgE',
            verbose=True,
        )

        assert data.check_completeness(verbose=False)['tracking'], 'Data object is missing attributes necessary for tracking.'
        return data

    def fetch_fov_list(self, sessions_to_exclude=None):
        """Fetch the list of FOVs to process."""
        match self.location:
            case 'remote':
                assert self.one and not self.one.offline, 'Remote data handling requires an active One instance'
                eids = self.one.search(subject=self.subject_path.name, datasets=['_ibl_mpciMeanImage.mlapdv_estimate.npy'])
                raise NotImplementedError('Remote data handling not implemented yet')
        # Find ordered list of 
        success, files, _ = self.input_files[0].find_files(self.subject_path)
        # Filter sessions with both meanimage and raw zip present
        files_grouped = groupby(sorted(files, key=lambda x: x.parent), key=lambda x: x.parent)
        alf_paths = [ALFPath(k) for k, g in files_grouped if len(list(g)) == 2]  # both files present
        
        if sessions_to_exclude:
            alf_paths = [path for path in alf_paths if path.session_path_short() not in sessions_to_exclude]
        return alf_paths
    
    def group_fovs(self, alf_paths):
        """Group FOVs by their location.
        
        NB: This can be simplified when we've accounted for objective angle
        """
    
        def most_common_suffix(paths):
            """
            Finds the most common numerical string suffix from a list of file paths.
            Assumes each path ends with a numerical string (e.g., '00', '04').
            
            Parameters:
            paths (list of str): List of file paths.
            
            Returns:
            str: Most common numerical suffix.
            """
            # Extract the suffix from each path (last 2 characters assumed)
            suffixes = [path.name[-2:] for path in paths]
            
            # Count occurrences of each suffix
            count = Counter(suffixes)
            
            # Find the most common one
            most_common = count.most_common(1)
            
            return most_common[0][0] if most_common else None

        # Exclude paths without mlapdv estimates
        alf_paths = sorted(p for p in alf_paths if (p / 'mpciMeanImage.mlapdv_estimate.npy').exists())

        # figure out the candidate repeated sites using the MLAPDV data
        MLAPDVs = np.array([np.load(p / 'mpciMeanImage.mlapdv_estimate.npy') for p in alf_paths]) # shape (n, 512, 512, 3)  [coords in microns: (ML, AP, DV)]
        MLAPDVs_means = np.mean(MLAPDVs, axis=(1,2))

        # Choose your strict threshold in microns:
        THRESH = 300.0  # use <= THRESH

        # Complete-linkage clustering enforces max intra-cluster distance ≤ THRESH
        Z = linkage(MLAPDVs_means, method='complete', metric='euclidean')
        cLabels = fcluster(Z, t=THRESH, criterion='distance')  # 1..K

        # Turn singletons into noise (-1)
        unique, counts = np.unique(cLabels, return_counts=True)
        sizes = dict(zip(unique, counts))
        cLabels = np.array([lab if sizes[lab] >= 2 else -1 for lab in cLabels])

        cluster_ids = [lab for lab in np.unique(cLabels) if lab != -1]

        # --- Order clusters by size (desc) and REMAP cLabels to 0..N-1 ---
        cluster_ids, cluster_sizes = np.unique(cLabels[cLabels != -1], return_counts=True)
        # cluster_ids are the OLD ids; sort them by size descending
        order = np.argsort(-cluster_sizes)
        sorted_old_ids = cluster_ids[order]

        # Map old id -> new id (size rank). Largest becomes 0.
        old2new = {old: new for new, old in enumerate(sorted_old_ids)}
        cLabels = np.array([old2new.get(lab, -1) for lab in cLabels])

        # Recompute sizes on the NEW cLabels (0..N-1)
        new_ids, new_sizes = np.unique(cLabels[cLabels != -1], return_counts=True)
        n_clusters = len(new_ids)

        clusters = dict.fromkeys(map(int, new_ids))

        for cID in new_ids:
            #in cases where two clustered FOVs come from the same session, choose the top one (TODO this is a HACK for now)
            path_list = []
            for p, n in zip(alf_paths, cLabels):
                if n == cID and (not any(path_list) or p.session_path() != path_list[-1].session_path()):
                    path_list.append(p)

            clusters[int(cID)] = path_list
            if path_list:
                #in cases where two clustered FOVs come from the same session, choose the top one (TODO this is a HACK for now)
                inferred_name = 'FOV_'+most_common_suffix([str(p)[-25:] for p in np.array(FOV_paths_sel)[idxs]])
                cluster_name = f"cFOV_{cID:02d}"
                
                logger.info('\nCluster name: %s   Inferred name: %s', cluster_name, inferred_name)
                [logger.info(f"{i:02d} {path}") for i,path in enumerate(path_list)]

                #[print(path[34:-11]) for i,path in enumerate(paths_toFOVs)]

        return clusters



if __name__ == '__main__':
    from pathlib import Path
    ROOT = Path('/mnt/whiterussian/Subjects')
    subject = 'SP072'
    subject_path = ROOT / subject
    
    task = ROICaTTask(subject_path)
    task.get_signatures()
    task.assert_expected_inputs()  # TODO support subject path
    
    # FOV_paths_all = list(np.unique([(Path(path) / '..').resolve() for path in paths_allMLAPDV]))  # parent
    # RFM_paths_all = list(np.unique([(Path(path) / '..').resolve() for path in paths_allRFM]))

    # if select_RFM_days_only:
    #     alf_paths_1 = list(np.unique([(Path(path) / '..').resolve() for path in FOV_paths_all]))
    #     alf_paths_2 = list(np.unique([(Path(path) / '..').resolve() for path in RFM_paths_all]))
    #     alf_paths_full = list(np.unique(list(set(alf_paths_1).intersection(set(alf_paths_2)))))
    # else:
    #     alf_paths_full = list(np.unique([(Path(path) / '..').resolve() for path in FOV_paths_all]))

    # #exclude the buggy days
    # alf_paths = list([Path(path).resolve() for path in alf_paths_full if not any(d in path.__str__() for d in days_to_exclude)])
    exclude = {'SP072/2025-09-17/001', 'SP072/2025-09-18/001', 'SP072/2025-09-19/001'}  # These MLAPDV files have odd dimensions; (512, 756, 3), (512, 752, 3)
    task._run(sessions_to_exclude=exclude)
