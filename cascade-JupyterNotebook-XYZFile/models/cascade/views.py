import os
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, Http404
import datetime
from .apply import predict_NMR, preprocess, RBFSequence

#NN model
import sys,io,os
import pandas as pd
import numpy as np
import gzip, pickle, argparse, warnings
import pickle
import math

from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ForwardSDMolSupplier
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem

from itertools import islice

from nfp.preprocessing import MolAPreprocessor, GraphSequence
from .genConf import genConf

import keras
import keras.backend as K

from keras.callbacks import ModelCheckpoint, CSVLogger, LearningRateScheduler

from keras.layers import (Input, Embedding, Dense, BatchNormalization,
                                 Concatenate, Multiply, Add)

from keras.models import Model, load_model

from nfp.layers import (MessageLayer, GRUStep, Squeeze, EdgeNetwork,
                               ReduceBondToPro, ReduceBondToAtom,
                               GatherAtomToBond, ReduceAtomToPro)
from nfp.models import GraphModel
from django.views.decorators.csrf import csrf_exempt

from .valid import validate_smiles

from .news import news

import json
import redis
import uuid
import ast
#Load NN model

def home(request):
    context = {
        'news': news
    }
    return render(request, 'cascade/home.html', context)

def about(request):
    return render(request, 'cascade/about.html', {'title': 'About'})

def predict(request):
    return render(request, 'cascade/predict.html')

@csrf_exempt
def JSmol_startup(request):
    m = startup_mols[0]
    c = m.GetConformers()[0]

    coords = ''
    for i,a in enumerate(m.GetAtoms()):
        ixyz = c.GetAtomPosition(i)
        coords += "{} {} {} {}|".format(atoms[a.GetAtomicNum()], *ixyz)

    jsmol_command = "data \"model example\"|{}|testing|{}end \"model example\";show data".format(m.GetNumAtoms(), coords)

    return HttpResponse(jsmol_command)

@csrf_exempt
def Shift_startup(request):
    responseTxt = ''
    for _,r in startup_shift.iterrows():
        responseTxt += '%s,%s;' % (int(r['atom_index']), r['Shift'])

    return HttpResponse(responseTxt)

redis_client = redis.StrictRedis(host="localhost", port=6379, db=0, decode_responses=True)

@csrf_exempt
def predict_NMR(request):
    #get smile string from the client

    smiles = request.POST["smiles"]

    if not validate_smiles(smiles):
        message = "Input molecule is not allowed"
        return JsonResponse({'message':message, 'task_id':None})

    #create task id
    task_id = uuid.uuid4().hex
    detail_key = "task_detail_{}".format(task_id)
    #send smiles to redis client
    redis_client.set(detail_key, json.dumps({
        "smiles": smiles
    }))
    #push task to the queue
    redis_client.rpush("task_queue", task_id)

    message = "Molecule has been submitted to the queue"
    return JsonResponse({'message':message, 'task_id':task_id})


@csrf_exempt
def check_task(request):
    task_id = request.GET["task_id"]
    result_key = "task_result_{}".format(task_id)
    result = redis_client.get(result_key)

    if result:
        result = ast.literal_eval(result)
        #draw molecules
        drawer = rdMolDraw2D.MolDraw2DSVG(1000, 600)
        drawer.SetFontSize(.6)

        opts = drawer.drawOptions()

        weightedShift = {int(item.split(',')[0]): item.split(',')[1] for item in filter(None, result['weightedShiftTxt'].split(';'))}
        for k,v in weightedShift.items():
            opts.atomLabels[k-1] = v
        opts.clearBackground=False
        opts.bondLineWidth=0.5
        opts.padding=0.1
        opts.additionalAtomLabelPadding=0.3

        smiles = result['smiles']
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '').replace(':svg', '')

        result_html = render(request, "cascade/results.html",
                                {'jsmol_command': result['jsmol_command'],
                                 'svg': svg,
                                 'smiles': smiles,
                                 'weightedShift': result['weightedShiftTxt'],
                                 'confShift': result['confShiftTxt'],
                                 'relative_E': result['relative_E'],
                                 'taskId': task_id,})
        return result_html
    else:
        return HttpResponse('running')

def download(request, taskId):
    file_path = os.path.join("cascade", "results", "{}.tar.gz".format(taskId))
    if os.path.exists(file_path):
        with open(file_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="application/tar+gzip")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(file_path)
            return response
    raise Http404
