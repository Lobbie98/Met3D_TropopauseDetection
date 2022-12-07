/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2018 Marc Rautenhaus
**  Copyright 2018 Michael Kern
**
**  Computer Graphics and Visualization Group
**  Technische Universitaet Muenchen, Garching, Germany
**
**  Met.3D is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Met.3D is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Met.3D.  If not, see <http://www.gnu.org/licenses/>.
**
*******************************************************************************/

#include "indexbuffer.h"

#include "util/mutil.h"
#include "util/mexception.h"
#include "gxfw/mglresourcesmanager.h"

using namespace GL;

/******************************************************************************
***                     CONSTRUCTOR / DESTRUCTOR                            ***
*******************************************************************************/

MIndexBuffer::MIndexBuffer(Met3D::MDataRequest requestKey,
                           const uint32_t numElements)
        : MAbstractGPUDataItem(requestKey),
          indexBufferObject(0), numIndices(numElements)

{
}


MIndexBuffer::~MIndexBuffer()
{
    // Delete the vertex buffer object. If the VBO has never been created,
    // "vertexBufferObject" takes the value 0 (constructor!)); glDeleteBuffers
    // ignores 0 buffers.
    glDeleteBuffers(1, &indexBufferObject); CHECK_GL_ERROR;
}


/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/


unsigned int MIndexBuffer::getGPUMemorySize_kb()
{
    return numIndices * sizeof(uint32_t) / 1024;
}


void MIndexBuffer::bindToElementArrayBuffer()
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBufferObject);
}


void MIndexBuffer::upload(const QVector<uint32_t>& data, QGLWidget* currentGLContext)
{
    upload(data.data(), static_cast<uint32_t>(data.size()), currentGLContext);
}


void MIndexBuffer::upload(const GLuint* data, const uint32_t elemCount,
                          QGLWidget* currentGLContext)
{
    // Make sure that the GLResourcesManager context is the currently active
    // context, otherwise glDrawArrays on the IBO generated here will fail in
    // any other context than the currently active. The GLResourcesManager
    // context is shared with all visible contexts.
    Met3D::MGLResourcesManager *glRM = Met3D::MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    // Delete the old IBO. If this is the first time the method has been called
    // (i.e. no IBO exists), "ibo" is set to 0 (constructor!); glDeleteBuffers
    // ignores 0 buffers.
    glDeleteBuffers(1, &indexBufferObject); CHECK_GL_ERROR;

    // Generate new VBO and upload geometry data to GPU.
    glGenBuffers(1, &indexBufferObject); CHECK_GL_ERROR;

    LOG4CPLUS_DEBUG(mlog, "size of GLuint: " << sizeof(GLuint) << " | size of uint32_t: " << sizeof(uint32_t));

    LOG4CPLUS_TRACE(mlog, "uploading index buffer data to ibo #"
            << indexBufferObject << flush);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBufferObject); CHECK_GL_ERROR;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 sizeof(uint32_t) * elemCount,
                 data,
                 GL_STATIC_DRAW); CHECK_GL_ERROR;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); CHECK_GL_ERROR;

    // If a valid GL context has been specifed, make this current.
    if (currentGLContext) currentGLContext->makeCurrent();
}
