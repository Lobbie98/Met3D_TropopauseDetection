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

#ifndef MET_3D_INDEXBUFFER_H
#define MET_3D_INDEXBUFFER_H

#include <log4cplus/loggingmacros.h>

// related third party imports
#include "GL/glew.h"
#include "QtCore"
#include "QGLWidget"


// local application imports
#include "gxfw/gl/abstractgpudataitem.h"

namespace GL
{

class MIndexBuffer : public MAbstractGPUDataItem
{
public:
    explicit MIndexBuffer(Met3D::MDataRequest requestKey, const uint32_t numElements);
    virtual ~MIndexBuffer();

    unsigned int getGPUMemorySize_kb() override;

    GLuint getIndexBufferObject() { return indexBufferObject; }

    void bindToElementArrayBuffer();

    void upload(const QVector<GLuint>& data, QGLWidget* currentGLContext = 0);
    void upload(const GLuint* data, const uint32_t elemCount,
                QGLWidget* currentGLContext = 0);

protected:
    GLuint indexBufferObject;
    GLuint numIndices;
};

}


#endif //MET_3D_INDEXBUFFER_H
