import numpy as np
import os


def read_skeleton(file):
    with open(file, 'r') as f:
        skeleton_sequence = {}
        skeleton_sequence['numFrame'] = int(f.readline())
        skeleton_sequence['frameInfo'] = []
        for t in range(skeleton_sequence['numFrame']):
            frame_info = {}
            frame_info['numBody'] = int(f.readline())
            frame_info['bodyInfo'] = []
            for m in range(frame_info['numBody']):
                body_info = {}
                body_info_key = [
                    'bodyID', 'clipedEdges', 'handLeftConfidence',
                    'handLeftState', 'handRightConfidence', 'handRightState',
                    'isResticted', 'leanX', 'leanY', 'trackingState'
                ]
                body_info = {
                    k: float(v)
                    for k, v in zip(body_info_key,
                                    f.readline().split())
                }
                body_info['numJoint'] = int(f.readline())
                body_info['jointInfo'] = []
                for v in range(body_info['numJoint']):
                    joint_info_key = [
                        'x', 'y', 'z', 'depthX', 'depthY', 'colorX', 'colorY',
                        'orientationW', 'orientationX', 'orientationY',
                        'orientationZ', 'trackingState'
                    ]
                    joint_info = {
                        k: float(v)
                        for k, v in zip(joint_info_key,
                                        f.readline().split())
                    }
                    body_info['jointInfo'].append(joint_info)
                frame_info['bodyInfo'].append(body_info)
            skeleton_sequence['frameInfo'].append(frame_info)
    return skeleton_sequence


# added by Tang 11-10
def read_rgb(filergb):
    with open(filergb, 'r') as ff:
        rgb_sequence = {}
        rgb_sequence['numFrame'] = int(ff.readline())
        rgb_sequence['frameInfo'] = []
        for t in range(rgb_sequence['numFrame']):
            frame_info = {}
            frame_info['numBody'] = int(ff.readline())
            frame_info['bodyInfo'] = []
            for m in range(frame_info['numBody']):
                body_info = {}
                body_info['numJoint'] = int(ff.readline())
                body_info['jointInfo'] = []
                for v in range(body_info['numJoint']):
                    joint_info_key = [
                        'r', 'g', 'b'
                    ]
                    joint_info = {
                        k: float(v)
                        for k, v in zip(joint_info_key,
                                        ff.readline().split())
                    }
                    body_info['jointInfo'].append(joint_info)
                frame_info['bodyInfo'].append(body_info)
            rgb_sequence['frameInfo'].append(frame_info)
    return rgb_sequence


'''
def read_xyz(file, max_body=2, num_joint=25):
    seq_info = read_skeleton(file)
    data = np.zeros((3, seq_info['numFrame'], num_joint, max_body))
    for n, f in enumerate(seq_info['frameInfo']):
        for m, b in enumerate(f['bodyInfo']):
            for j, v in enumerate(b['jointInfo']):
                if m < max_body and j < num_joint:
                    data[:, n, j, m] = [v['x'], v['y'], v['z']]
                else:
                    pass
    return data
''' 

# added by Tang 11-10
def read_xyz(file, filergb, max_body=2, num_joint=25):
    '''
    seq_info = read_skeleton(file)
    rgb_info = read_rgb(filergb)

    data = np.zeros((6, seq_info['numFrame'], num_joint, max_body))
    for n, f in enumerate(seq_info['frameInfo']):
        for m, b in enumerate(f['bodyInfo']):
            for j, v in enumerate(b['jointInfo']):
                if m < max_body and j < num_joint:
                    data[0, n, j, m] = v['x']
                    data[1, n, j, m] = v['y']
                    data[2, n, j, m] = v['z']
                else:
                    pass

    for n, f in enumerate(rgb_info['frameInfo']):
        nn=(n-1)*10+1
        for m, b in enumerate(f['bodyInfo']):
            for j, v in enumerate(b['jointInfo']):
                if m < max_body and j < num_joint:

                    data[3, nn, j, m] = v['r']
                    data[4, nn, j, m] = v['g']
                    data[5, nn, j, m] = v['b']

                    nextn=n*10+1
                    if nextn>seq_info['numFrame']:
                        kk=nn
                        while(kk<seq_info['numFrame']):
                            data[3, kk, j, m] = v['r']
                            data[4, kk, j, m] = v['g']
                            data[5, kk, j, m] = v['b']
                            kk=kk+1
                    else:
                        kk=1
                        while(kk<10):
                            data[3, nn+kk, j, m] = v['r']
                            data[4, nn+kk, j, m] = v['g']
                            data[5, nn+kk, j, m] = v['b']
                            kk=kk+1
                else:
                    pass
    '''

    seq_info = read_skeleton(file)
    #rgb_info = read_rgb(filergb)
    if os.path.exists(filergb):
        rgb_info = read_rgb(filergb)
        data = np.zeros((6, seq_info['numFrame'], num_joint, max_body))
        for n, f in enumerate(seq_info['frameInfo']):
            for m, b in enumerate(f['bodyInfo']):
                for j, v in enumerate(b['jointInfo']):
                    if m < max_body and j < num_joint:
                        data[0, n, j, m] = v['x']
                        data[1, n, j, m] = v['y']
                        data[2, n, j, m] = v['z']
                    else:
                        pass

        for n, f in enumerate(rgb_info['frameInfo']):
            nn=(n-1)*10+1
            for m, b in enumerate(f['bodyInfo']):
                for j, v in enumerate(b['jointInfo']):
                    if m < max_body and j < num_joint:

                        data[3, nn, j, m] = v['r']
                        data[4, nn, j, m] = v['g']
                        data[5, nn, j, m] = v['b']

                        nextn=n*10+1
                        if nextn>seq_info['numFrame']:
                            kk=nn
                            while(kk<seq_info['numFrame']):
                                data[3, kk, j, m] = v['r']
                                data[4, kk, j, m] = v['g']
                                data[5, kk, j, m] = v['b']
                                kk=kk+1
                        else:
                            kk=1
                            while(kk<10):
                                data[3, nn+kk, j, m] = v['r']
                                data[4, nn+kk, j, m] = v['g']
                                data[5, nn+kk, j, m] = v['b']
                                kk=kk+1
                    else:
                        pass
    else:
        data = np.zeros((6, seq_info['numFrame'], num_joint, max_body))
        for n, f in enumerate(seq_info['frameInfo']):
            for m, b in enumerate(f['bodyInfo']):
                for j, v in enumerate(b['jointInfo']):
                    if m < max_body and j < num_joint:
                        data[0, n, j, m] = v['x']
                        data[1, n, j, m] = v['y']
                        data[2, n, j, m] = v['z']
                        data[3, n, j, m] = 0
                        data[4, n, j, m] = 0
                        data[5, n, j, m] = 0
                    else:
                        pass

    return data
